const DOTPLOT_RENDERING = {
  intersectionMinSizePx: 2
};

// Normalized bounds of the real plotting area inside the dotplot SVG.
// Values are ratios of the displayed SVG box (0 = left/top, 1 = right/bottom).
// Adjust these by trial and error to align external genomic tracks with SVG axes.
const DOTPLOT_AXIS_BOUNDS = {
  // X ratios: fraction of SVG width, measured from the left (0 = left, 1 = right).
  // xZeroRatio = pixel position of genomic coordinate 1 (left end of x axis).
  // xMaxRatio  = pixel position of maximum genomic coordinate (right end of x axis).
  xZeroRatio: 0, //0.05,
  xMaxRatio: 1, //0.992,
  // Y ratios: measured from the BOTTOM of the SVG (0 = bottom, 1 = top).
  // yZeroRatio = position of genomic coordinate 1 (bottom of y axis).
  // yMaxRatio  = position of maximum genomic coordinate (top of y axis).
  // Conversion to CSS pixel y: pixelY = imageHeight * (1 - ratio).
  yZeroRatio: 0, //0.0668,
  yMaxRatio: 1 //0.9683
};

// Set to true to draw calibration lines at axis boundaries on track canvases.
const DOTPLOT_DEBUG_LAYOUT = false;

// Track dimensions include one feature inset on each side so the inner region
// matches the browser-mode sample track.
const DOTPLOT_TRACK = {
  yTrackWidth:    TRACK_GEOMETRY.trackHeight + 2 * TRACK_GEOMETRY.featureInset,
  xTrackHeight:   TRACK_GEOMETRY.trackHeight + 2 * TRACK_GEOMETRY.featureInset,
  debugColor:     "#3b82f6",
  debugLineWidth: 1.5
};

// TRACK_HIGHLIGHT_INSET controls highlight stroke positioning independently of
// TRACK_GEOMETRY.featureInset, which positions rendered feature geometry.
const TRACK_HIGHLIGHT_INSET = 0.5;

// Dotplot zoom settings.
const DOTPLOT_ZOOM_STEP = 1.3;
const DOTPLOT_ZOOM_MIN  = 0.25;
const DOTPLOT_ZOOM_MAX  = 10;

const DOTPLOT_AXIS = {
  tickLength: 5,
  fontSize: 10,
  labelPadding: 2,
  color: "#444444",
  labelColor: "#555555",
  targetTickSpacingPx: 90
};

// The shared matrix/X/Y genomic origin lives at local (0, 0), so the context
// axes meet the matrix boundary without a corner gap.  Only the right margin
// remains for X-axis label rendering near the image end.
const DOTPLOT_OUTER_PADDING = {
  top: 0,
  right: 10,
  bottom: 0,
  left: 0
};

// This is frame-only space between the outer dotplot border and the scrollport.
// It deliberately does not participate in any Konva coordinate system.
const DOTPLOT_FRAME_TOP_BREATHING_SPACE = 5;

// Geometry for the currently rendered dotplot stage. Projection overlays must
// share this object with the image and tracks rather than sampling DOM width
// again after the stage has changed the scroll layout.
let _currentDotplotGeometry = null;

// Pair changes may replace the SVG with one that has a slightly different
// fitted geometry. Keep a one-shot genomic viewport center so the new pair can
// restore that center after its three panes have been laid out.
let _pendingDotplotViewportRestore = null;
let _dotplotViewportRestoreEpoch = 0;
let _renderedDotplotPairKey = null;

function getDotplotPairs() {
  return (REGION_DATA.dotplots && REGION_DATA.dotplots.pairs) || [];
}

function buildDotplotPairIndexes() {
  dotplotPairIndexes.validPairsByY.clear();
  dotplotPairIndexes.validPairsByX.clear();
  dotplotPairIndexes.pairByYX.clear();

  const samples = new Set();

  for (const pair of getDotplotPairs()) {
    const ySample = pair.y_sample;
    const xSample = pair.x_sample;

    samples.add(ySample);
    samples.add(xSample);

    if (!dotplotPairIndexes.validPairsByY.has(ySample)) {
      dotplotPairIndexes.validPairsByY.set(ySample, new Set());
    }
    dotplotPairIndexes.validPairsByY.get(ySample).add(xSample);

    if (!dotplotPairIndexes.validPairsByX.has(xSample)) {
      dotplotPairIndexes.validPairsByX.set(xSample, new Set());
    }
    dotplotPairIndexes.validPairsByX.get(xSample).add(ySample);

    dotplotPairIndexes.pairByYX.set(`${ySample}::${xSample}`, pair);
  }

  dotplotPairIndexes.samples = [...samples].sort();
}

function findDotplotPair(ySample, xSample) {
  if (!ySample || !xSample) {
    return null;
  }

  return dotplotPairIndexes.pairByYX.get(`${ySample}::${xSample}`) || null;
}

function getSampleByName(sampleName) {
  return searchIndexes.sampleByName.get(sampleName) || null;
}

// Returns the fixed footprint of the GFF strips only.  The categorical legend
// is viewport DOM UI, deliberately outside the genomic X stage.
function getDotplotYGffTotalWidth(ySampleData) {
  const trackCount = ySampleData ? getSampleGffTracks(ySampleData).length : 0;
  return trackCount * (GFF_TRACK.height + GFF_TRACK.gap);
}

function getDotplotXGffTotalHeight(xSampleData) {
  const trackCount = xSampleData ? getSampleGffTracks(xSampleData).length : 0;
  return trackCount > 0
    ? GFF_TRACK.topGap + trackCount * (GFF_TRACK.height + GFF_TRACK.gap)
    : 0;
}

function getDotplotLayoutViewport() {
  const container = document.querySelector(".dotplot-content");
  const viewportWidth = Math.max(1, container ? container.clientWidth : 800);
  return { viewportWidth };
}

function getDotplotImageDisplaySize(img, maxWidth) {
  const naturalWidth = Math.max(1, img.naturalWidth);
  const naturalHeight = Math.max(1, img.naturalHeight);
  // The zoom=1 frame is width-responsive. Zoom is applied only after this
  // intrinsic-aspect-ratio base size has been established.
  const fitScale = maxWidth / naturalWidth;
  const baseImageWidth = Math.max(1, Math.round(naturalWidth * fitScale));
  const baseImageHeight = Math.max(1, Math.round(naturalHeight * fitScale));

  return {
    baseImageWidth,
    baseImageHeight,
    imageWidth: Math.max(1, Math.round(baseImageWidth * _dotplotState.zoom)),
    imageHeight: Math.max(1, Math.round(baseImageHeight * _dotplotState.zoom))
  };
}

// Computes every local coordinate system from one logical layout calculation.
// Matrix and X context share the same X genomic pixels; matrix and Y context
// share the same Y genomic pixels.  Only the stage-local origins differ.
function computeDotplotGeometry(horizontalScrollbarHeight = 0) {
  const img = document.getElementById("dotplot-svg-img");
  if (!img || !img.complete || img.naturalWidth === 0) {
    return null;
  }

  const xSampleData = getSampleByName(_dotplotState.selectedX);
  const ySampleData = getSampleByName(_dotplotState.selectedY);
  const viewport = getDotplotLayoutViewport();
  const yGffWidth = getDotplotYGffTotalWidth(ySampleData);
  const yGffSideGap = yGffWidth > 0 ? GFF_TRACK.topGap : 0;
  const xGffHeight = getDotplotXGffTotalHeight(xSampleData);
  const xAxisGap = getDotplotXAxisGap();

  // Estimate the Y label gutter before fitting; its displayed units are based
  // on the fixed genomic span, so the result is stable for this redraw.
  const yAxisGap = getDotplotYAxisGap(ySampleData, Math.max(1, img.naturalHeight));
  const yContextWidth = DOTPLOT_OUTER_PADDING.left
    + yGffWidth + yGffSideGap + DOTPLOT_TRACK.yTrackWidth + yAxisGap;
  const xContextHeight = xAxisGap + DOTPLOT_TRACK.xTrackHeight + xGffHeight;

  const availableImageWidth = Math.max(
    100,
    viewport.viewportWidth - yContextWidth - DOTPLOT_OUTER_PADDING.left - DOTPLOT_OUTER_PADDING.right
  );
  const size = getDotplotImageDisplaySize(
    img,
    availableImageWidth
  );
  const { baseImageWidth, baseImageHeight, imageWidth, imageHeight } = size;
  const baseMatrixWidth = DOTPLOT_OUTER_PADDING.left + baseImageWidth + DOTPLOT_OUTER_PADDING.right;
  const baseMatrixHeight = DOTPLOT_OUTER_PADDING.top + baseImageHeight + DOTPLOT_OUTER_PADDING.bottom;
  const matrixWidth = DOTPLOT_OUTER_PADDING.left + imageWidth + DOTPLOT_OUTER_PADDING.right;
  const matrixHeight = DOTPLOT_OUTER_PADDING.top + imageHeight + DOTPLOT_OUTER_PADDING.bottom;
  const imageX = DOTPLOT_OUTER_PADDING.left;
  const imageY = DOTPLOT_OUTER_PADDING.top;
  const xZero = imageX + imageWidth * DOTPLOT_AXIS_BOUNDS.xZeroRatio;
  const xMax = imageX + imageWidth * DOTPLOT_AXIS_BOUNDS.xMaxRatio;
  const yZeroPixel = imageY + imageHeight * (1 - DOTPLOT_AXIS_BOUNDS.yZeroRatio);
  const yMaxPixel = imageY + imageHeight * (1 - DOTPLOT_AXIS_BOUNDS.yMaxRatio);
  // The zoom=1 fit fixes the frame footprint. When a horizontal scrollbar is
  // present, it consumes space inside that footprint rather than enlarging it.
  const baseMatrixViewportHeight = baseMatrixHeight;
  const scrollportHeight = baseMatrixViewportHeight + xContextHeight;
  const matrixViewportHeight = Math.min(
    matrixHeight,
    Math.max(1, baseMatrixViewportHeight - horizontalScrollbarHeight)
  );

  const xContext = {
    width: matrixWidth,
    height: xContextHeight,
    xZero,
    xMax,
    axisY: 0.5,
    trackBoxY: xAxisGap,
    outerPadding: DOTPLOT_OUTER_PADDING,
    xGffHeight
  };
  const yContext = {
    width: yContextWidth,
    height: matrixHeight,
    yZeroPixel,
    yMaxPixel,
    axisX: yContextWidth - 0.5,
    trackBoxX: DOTPLOT_OUTER_PADDING.left + yGffWidth + yGffSideGap,
    outerPadding: DOTPLOT_OUTER_PADDING,
    yGffWidth,
    yGffSideGap
  };

  return {
    baseImageWidth,
    baseImageHeight,
    baseMatrixWidth,
    baseMatrixHeight,
    matrixWidth,
    matrixHeight,
    imageX,
    imageY,
    imageWidth,
    imageHeight,
    xZero,
    xMax,
    yZeroPixel,
    yMaxPixel,
    outerPadding: DOTPLOT_OUTER_PADDING,
    xAxisGap,
    yAxisGap,
    xContext,
    yContext,
    matrixViewportWidth: Math.min(matrixWidth, Math.max(1, viewport.viewportWidth - yContextWidth)),
    matrixViewportHeight,
    horizontalScrollbarHeight,
    viewportWidth: viewport.viewportWidth,
    scrollportWidth: Math.min(viewport.viewportWidth, yContextWidth + matrixWidth),
    scrollportHeight,
    surfaceWidth: yContextWidth + matrixWidth,
    surfaceHeight: matrixHeight + xContextHeight
  };
}

// Maps a genomic region position to Konva stage x for the x-sample track.
// position 1 → xZero; position region_length → xMax.
function mapXCoordinateToStagePx(position, sample, geometry) {
  const regionLength = sample.region_length;
  if (regionLength <= 1) {
    return geometry.xZero;
  }
  const ratio = (position - 1) / (regionLength - 1);
  return geometry.xZero + ratio * (geometry.xMax - geometry.xZero);
}

// Maps a genomic region position to Konva stage y for the y-sample track.
// position 1 → yZeroPixel (near image bottom); position region_length → yMaxPixel (near image top).
function mapYCoordinateToStagePx(position, sample, geometry) {
  const regionLength = sample.region_length;
  if (regionLength <= 1) {
    return geometry.yZeroPixel;
  }
  const ratio = (position - 1) / (regionLength - 1);
  return geometry.yZeroPixel - ratio * (geometry.yZeroPixel - geometry.yMaxPixel);
}

function formatAxisValueForSpan(value, visibleSpan) {
  return formatGenomicCoordinate(value, visibleSpan);
}

function getDotplotAxisTickValues(sample, pixelSpan) {
  if (!sample) {
    return [];
  }

  const visibleStart = 1;
  const visibleEnd = sample.region_length;
  const visibleSpan = Math.max(1, visibleEnd - visibleStart + 1);
  const bpPerPixel = visibleSpan / Math.max(1, pixelSpan);
  const rawStep = bpPerPixel * DOTPLOT_AXIS.targetTickSpacingPx;
  const step = niceStep(rawStep);
  const firstTick = Math.ceil(visibleStart / step) * step;
  const values = [];

  for (let value = firstTick; value <= visibleEnd; value += step) {
    values.push(value);
  }

  return values;
}


function estimateDotplotAxisLabelWidth(sample, pixelSpan) {
  if (!sample) {
    return 0;
  }

  const visibleSpan = Math.max(1, sample.region_length);
  const tickValues = getDotplotAxisTickValues(sample, pixelSpan);

  if (tickValues.length === 0) {
    return 0;
  }

  return Math.max(
    ...tickValues.map(value =>
      estimateTextWidth(
        formatAxisValueForSpan(value, visibleSpan),
        DOTPLOT_AXIS.fontSize
      )
    )
  );
}


function getDotplotXAxisGap() {
  return Math.ceil(
    DOTPLOT_AXIS.tickLength
    + DOTPLOT_AXIS.labelPadding
    + DOTPLOT_AXIS.fontSize
    + 8
  );
}


function getDotplotYAxisGap(ySampleData, imageHeight) {
  const labelWidth = estimateDotplotAxisLabelWidth(ySampleData, imageHeight);

  return Math.ceil(
    labelWidth
    + DOTPLOT_AXIS.tickLength
    + DOTPLOT_AXIS.labelPadding
    + 10
  );
}

function drawDotplotXAxis(layer, geometry, sample) {
  if (!sample) {
    return;
  }

  const axisY = geometry.axisY;
  const visibleStart = 1;
  const visibleEnd = sample.region_length;
  const visibleSpan = Math.max(1, visibleEnd - visibleStart + 1);

  const bpPerPixel = visibleSpan / Math.max(1, geometry.xMax - geometry.xZero);
  const rawStep = bpPerPixel * DOTPLOT_AXIS.targetTickSpacingPx;
  const step = niceStep(rawStep);
  const firstTick = Math.ceil(visibleStart / step) * step;

  layer.add(new Konva.Line({
    points: [geometry.xZero, axisY, geometry.xMax, axisY],
    stroke: DOTPLOT_AXIS.color,
    strokeWidth: 1,
    listening: false
  }));

  for (let value = firstTick; value <= visibleEnd; value += step) {
    const x = mapXCoordinateToStagePx(value, sample, geometry);

    layer.add(new Konva.Line({
      points: [x, axisY, x, axisY + DOTPLOT_AXIS.tickLength],
      stroke: DOTPLOT_AXIS.color,
      strokeWidth: 1,
      listening: false
    }));

    layer.add(new Konva.Text({
      x: x - 34,
      y: axisY + DOTPLOT_AXIS.tickLength + DOTPLOT_AXIS.labelPadding,
      width: 68,
      text: formatAxisValueForSpan(value, visibleSpan),
      fontSize: DOTPLOT_AXIS.fontSize,
      fill: DOTPLOT_AXIS.labelColor,
      align: "center",
      listening: false
    }));
  }
}

function drawDotplotYAxis(layer, geometry, sample) {
  if (!sample) {
    return;
  }

  const axisX = geometry.axisX;
  const visibleStart = 1;
  const visibleEnd = sample.region_length;
  const visibleSpan = Math.max(1, visibleEnd - visibleStart + 1);

  const axisHeight = Math.max(1, geometry.yZeroPixel - geometry.yMaxPixel);
  const bpPerPixel = visibleSpan / axisHeight;
  const rawStep = bpPerPixel * DOTPLOT_AXIS.targetTickSpacingPx;
  const step = niceStep(rawStep);
  const firstTick = Math.ceil(visibleStart / step) * step;

  layer.add(new Konva.Line({
    points: [axisX, geometry.yZeroPixel, axisX, geometry.yMaxPixel],
    stroke: DOTPLOT_AXIS.color,
    strokeWidth: 1,
    listening: false
  }));

  for (let value = firstTick; value <= visibleEnd; value += step) {
    const y = mapYCoordinateToStagePx(value, sample, geometry);

    layer.add(new Konva.Line({
      points: [axisX - DOTPLOT_AXIS.tickLength, y, axisX, y],
      stroke: DOTPLOT_AXIS.color,
      strokeWidth: 1,
      listening: false
    }));

    layer.add(new Konva.Text({
      x: axisX - DOTPLOT_AXIS.tickLength - 70,
      y: y - DOTPLOT_AXIS.fontSize / 2,
      width: 66,
      text: formatAxisValueForSpan(value, visibleSpan),
      fontSize: DOTPLOT_AXIS.fontSize,
      fill: DOTPLOT_AXIS.labelColor,
      align: "right",
      listening: false
    }));
  }
}

// Returns true when dotplot mode is the active viewer mode.
function isDotplotModeActive() {
  const panel = document.getElementById("dotplot-panel");
  return panel ? !panel.classList.contains("hidden") : false;
}

function getSelectedDotplotPairKey() {
  const pair = findDotplotPair(_dotplotState.selectedY, _dotplotState.selectedX);
  return pair ? `${pair.y_sample}::${pair.x_sample}` : null;
}

function captureDotplotViewportCenter() {
  const geometry = _currentDotplotGeometry;
  const container = document.querySelector(".dotplot-content");
  if (!geometry || !container) {
    return null;
  }

  const xSpan = geometry.xMax - geometry.xZero;
  const ySpan = geometry.yZeroPixel - geometry.yMaxPixel;
  if (xSpan <= 0 || ySpan <= 0 || geometry.matrixViewportWidth <= 0 || geometry.matrixViewportHeight <= 0) {
    return null;
  }

  const matrixCenterX = container.scrollLeft + geometry.matrixViewportWidth / 2;
  const matrixCenterY = container.scrollTop + geometry.matrixViewportHeight / 2;
  return {
    xRatio: clampValue((matrixCenterX - geometry.xZero) / xSpan, 0, 1),
    yRatio: clampValue((geometry.yZeroPixel - matrixCenterY) / ySpan, 0, 1)
  };
}

function restorePendingDotplotViewportCenter() {
  const restore = _pendingDotplotViewportRestore;
  if (!restore) {
    return;
  }
  _pendingDotplotViewportRestore = null;

  // Let the browser commit the new canvas sizes before applying the one-shot
  // scroll position. This is deliberately not tied to native scroll events.
  requestAnimationFrame(() => {
    if (restore.epoch !== _dotplotViewportRestoreEpoch || restore.pairKey !== getSelectedDotplotPairKey()) {
      return;
    }

    const geometry = _currentDotplotGeometry;
    const container = document.querySelector(".dotplot-content");
    if (!geometry || !container) {
      return;
    }

    const matrixCenterX = geometry.xZero + restore.xRatio * (geometry.xMax - geometry.xZero);
    const matrixCenterY = geometry.yZeroPixel - restore.yRatio * (geometry.yZeroPixel - geometry.yMaxPixel);
    container.scrollLeft = clampValue(
      matrixCenterX - geometry.matrixViewportWidth / 2,
      0,
      Math.max(0, container.scrollWidth - container.clientWidth)
    );
    container.scrollTop = clampValue(
      matrixCenterY - geometry.matrixViewportHeight / 2,
      0,
      Math.max(0, container.scrollHeight - container.clientHeight)
    );
  });
}

function preserveDotplotViewportForNextRedraw() {
  if (!isDotplotModeActive()) {
    return;
  }

  const pairKey = getSelectedDotplotPairKey();
  const previousViewportCenter = pairKey ? captureDotplotViewportCenter() : null;
  if (!previousViewportCenter) {
    return;
  }

  const epoch = ++_dotplotViewportRestoreEpoch;
  _pendingDotplotViewportRestore = { ...previousViewportCenter, pairKey, epoch };
}

// Schedules a dotplot redraw via requestAnimationFrame.
// Guarded against duplicate pending frames and against firing in browser mode.
function requestDotplotRedraw() {
  if (_dotplotRedrawPending || !isDotplotModeActive()) {
    return;
  }
  _dotplotRedrawPending = true;
  requestAnimationFrame(() => {
    _dotplotRedrawPending = false;
    redrawDotplotStage();
  });
}

function createDotplotBlockHighlightShape(axis) {
  return new Konva.Shape({
    sceneFunc(ctx) {
      const geoms = axis === "x" ? _dotplotXBlockHighlightGeoms : _dotplotYBlockHighlightGeoms;
      const color = axis === "x" ? _dotplotXBlockHighlightColor : _dotplotYBlockHighlightColor;
      if (geoms.length === 0) { return; }
      ctx.save();
      ctx.fillStyle = color;
      ctx.beginPath();
      for (const rect of geoms) {
        ctx.rect(rect.x, rect.y, rect.width, rect.height);
      }
      ctx.fill();
      ctx.restore();
    },
    visible: false,
    listening: false
  });
}

function createDotplotSnpHighlightShape(axis) {
  return new Konva.Shape({
    sceneFunc(ctx) {
      const geoms = axis === "x" ? _dotplotXSnpHighlightGeoms : _dotplotYSnpHighlightGeoms;
      const color = axis === "x" ? _dotplotXSnpHighlightColor : _dotplotYSnpHighlightColor;
      if (geoms.length === 0) { return; }
      ctx.save();
      ctx.strokeStyle = color;
      ctx.lineWidth = FEATURE_RENDERING.snpHighlightMinWidthPx;
      ctx.beginPath();
      for (const line of geoms) {
        if (axis === "x") {
          ctx.moveTo(line.cx, line.y0);
          ctx.lineTo(line.cx, line.y1);
        } else {
          ctx.moveTo(line.x0, line.cy);
          ctx.lineTo(line.x1, line.cy);
        }
      }
      ctx.stroke();
      ctx.restore();
    },
    visible: false,
    listening: false
  });
}

function wireDotplotTrackInteraction(stage, axis) {
  stage.on("pointermove", () => {
    const pointer = stage.getPointerPosition();
    if (!pointer || state.isApplyingPin) { return; }
    const resolved = resolveDotplotHoveredFeature(axis, pointer.x, pointer.y);
    applyDotplotResolvedHover(resolved);

    const gene = resolved ? null : resolveDotplotGffGene(axis, pointer.x, pointer.y);
    if (gene) {
      showDotplotGffGeneTooltip(gene, stage, pointer.x, pointer.y);
    } else {
      hideGffGeneTooltip();
    }
    stage.container().style.cursor = resolved || gene ? "pointer" : "default";
  });

  stage.on("pointerleave", () => {
    stage.container().style.cursor = "default";
    hideGffGeneTooltip();
    if (!state.isApplyingPin) {
      applyDotplotResolvedHover(null);
    }
  });

  stage.on("click", () => {
    if (state.isApplyingPin) { return; }
    const pointer = stage.getPointerPosition();
    if (!pointer) { return; }
    const resolved = resolveDotplotHoveredFeature(axis, pointer.x, pointer.y);
    if (!resolved) { return; }
    state.isApplyingPin = true;
    state.hoveredFeatureType = null;
    state.hoveredFeatureId = null;
    _lastResolvedDotplotHoverKey = null;
    _lastResolvedHoverKey = null;
    setPinnedFeature(resolved.featureType, resolved.featureId);
    requestAnimationFrame(() => {
      state.isApplyingPin = false;
    });
  });
}

// Creates the three deliberately separate Dotplot surfaces.  Only the matrix
// stage contains the SVG; X and Y are lightweight annotation strips.
function initDotplotStage() {
  if (dotplotStage) { return; }

  dotplotStage = new Konva.Stage({ container: "dotplot-viewer", width: 1, height: 1 });
  dotplotImageLayer = new Konva.Layer({ listening: false });
  dotplotDebugLayer = new Konva.Layer({ listening: false });
  dotplotHighlightLayer = new Konva.Layer({ listening: false });
  dotplotStage.add(dotplotImageLayer, dotplotDebugLayer, dotplotHighlightLayer);

  dotplotXStage = new Konva.Stage({ container: "dotplot-x-viewer", width: 1, height: 1 });
  dotplotXTrackLayer = new Konva.Layer({ listening: false });
  dotplotXHighlightLayer = new Konva.Layer({ listening: false });
  dotplotXStage.add(dotplotXTrackLayer, dotplotXHighlightLayer);

  dotplotYStage = new Konva.Stage({ container: "dotplot-y-viewer", width: 1, height: 1 });
  dotplotYTrackLayer = new Konva.Layer({ listening: false });
  dotplotYHighlightLayer = new Konva.Layer({ listening: false });
  dotplotYStage.add(dotplotYTrackLayer, dotplotYHighlightLayer);

  _dotplotXBlockHighlightShape = createDotplotBlockHighlightShape("x");
  _dotplotXSnpHighlightShape = createDotplotSnpHighlightShape("x");
  _dotplotYBlockHighlightShape = createDotplotBlockHighlightShape("y");
  _dotplotYSnpHighlightShape = createDotplotSnpHighlightShape("y");

  _dotplotBlockIntersectionShape = new Konva.Shape({
    sceneFunc(ctx) {
      if (!_dotplotBlockIntersectionGeom) { return; }
      const { vertical, horizontal } = _dotplotBlockIntersectionGeom;
      ctx.save();
      ctx.fillStyle = "rgba(59, 130, 246, 0.18)";
      ctx.fillRect(vertical.x, vertical.y, vertical.width, vertical.height);
      ctx.fillRect(horizontal.x, horizontal.y, horizontal.width, horizontal.height);
      ctx.restore();
    },
    visible: false,
    listening: false
  });
  _dotplotSnpProjectionShape = new Konva.Shape({
    sceneFunc(ctx) {
      if (!_dotplotSnpProjectionGeom) { return; }
      const { x, y, xZero, yZeroPixel } = _dotplotSnpProjectionGeom;
      ctx.save();
      ctx.strokeStyle = "rgba(59, 130, 246, 0.45)";
      ctx.lineWidth = FEATURE_RENDERING.snpHighlightMinWidthPx;
      ctx.beginPath();
      ctx.moveTo(x, y); ctx.lineTo(x, yZeroPixel);
      ctx.moveTo(xZero, y); ctx.lineTo(x, y);
      ctx.stroke();
      ctx.restore();
    },
    visible: false,
    listening: false
  });

  wireDotplotTrackInteraction(dotplotXStage, "x");
  wireDotplotTrackInteraction(dotplotYStage, "y");
}

// Full batched redraw of the dotplot Konva stage.
// Uses one Konva.Shape per feature group — no one-node-per-feature.
// Visual style matches browser mode: white region, gray blocks, PASS/rejected SNPs,
// black outline.
// Also rebuilds the hover spatial index so hit-testing is always in sync with the layout.
// Computes along-axis geometry for all blocks and SNPs of one sample track,
// in a "local horizontal" coordinate system where the primary axis runs along
// the track and the cross-axis is the track height (TRACK_GEOMETRY.trackHeight).
//
// Returns:
//   fillRects:    [{along0, len, featureId}]  — block fill segments along the primary axis
//   snpPositions: [{along, featureId}]        — SNP pixel positions along the primary axis
//
// `mapper` maps a region position (genomic coordinate) to a stage pixel along
// the track's primary axis.  For the x-track, mapper = mapXCoordinateToStagePx.
// For the y-track, mapper = mapYCoordinateToStagePx (y-axis is inverted).
function buildTrackAlongAxisGeoms(blocks, snps, mapper) {
  const fillRects = [];
  for (const block of blocks) {
    const px0    = mapper(block.block_start_in_region);
    const px1    = mapper(block.block_end_in_region);
    const along0 = Math.min(px0, px1);
    const len    = Math.max(FEATURE_RENDERING.blockMinWidthPx, Math.abs(px1 - px0));
    fillRects.push({ along0, len, featureId: block.feature_id });
  }
  const snpPositions = [];
  for (const snp of snps) {
    if (!shouldDisplaySnp(snp.feature_id)) {
      continue;
    }
    snpPositions.push({ along: mapper(snp.pos_in_region), featureId: snp.feature_id });
  }
  return { fillRects, snpPositions };
}

function splitDotplotSnpEntriesByStatus(entries) {
  const pass = [];
  const rejected = [];

  for (const entry of entries) {
    if (isRejectedSnp(entry.featureId)) {
      rejected.push(entry);
    } else {
      pass.push(entry);
    }
  }

  return { pass, rejected };
}

function drawDotplotSnpLines(layer, entries, color, drawLine) {
  if (entries.length === 0) {
    return;
  }

  layer.add(new Konva.Shape({
    sceneFunc(ctx, shape) {
      ctx.beginPath();
      for (const entry of entries) {
        drawLine(ctx, entry);
      }
      ctx.fillStrokeShape(shape);
    },
    stroke: color,
    strokeWidth: FEATURE_RENDERING.snpMinWidthPx,
    listening: false
  }));
}

function drawDotplotFeatureTrack(layer, axis, localGeometry, sample) {
  const isX = axis === "x";
  const mapper = isX
    ? position => mapXCoordinateToStagePx(position, sample, localGeometry)
    : position => mapYCoordinateToStagePx(position, sample, localGeometry);
  const regionX = isX
    ? localGeometry.xZero
    : localGeometry.trackBoxX + TRACK_GEOMETRY.featureInset;
  const regionY = isX
    ? localGeometry.trackBoxY + TRACK_GEOMETRY.featureInset
    : localGeometry.yMaxPixel;
  const regionW = isX
    ? Math.max(1, localGeometry.xMax - localGeometry.xZero)
    : TRACK_GEOMETRY.trackHeight;
  const regionH = isX
    ? TRACK_GEOMETRY.trackHeight
    : Math.max(1, localGeometry.yZeroPixel - localGeometry.yMaxPixel);
  const geoms = sample
    ? buildTrackAlongAxisGeoms(sample.blocks, sample.snps, mapper)
    : { fillRects: [], snpPositions: [] };
  const blocks = geoms.fillRects.map(rect => isX
    ? { x: rect.along0, y: regionY + TRACK_GEOMETRY.featureInset, width: rect.len, height: getFeatureHeight(), featureId: rect.featureId }
    : { x: regionX + TRACK_GEOMETRY.featureInset, y: rect.along0, width: getFeatureHeight(), height: rect.len, featureId: rect.featureId }
  );
  const snps = geoms.snpPositions.map(entry => isX
    ? { cx: entry.along, y0: regionY + TRACK_GEOMETRY.featureInset, y1: regionY + regionH - TRACK_GEOMETRY.featureInset, featureId: entry.featureId }
    : { cy: entry.along, x0: regionX + TRACK_GEOMETRY.featureInset, x1: regionX + regionW - TRACK_GEOMETRY.featureInset, featureId: entry.featureId }
  );

  layer.add(new Konva.Rect({ x: regionX, y: regionY, width: regionW, height: regionH, fill: "#ffffff", listening: false }));
  if (blocks.length > 0) {
    layer.add(new Konva.Shape({
      sceneFunc(ctx, shape) {
        ctx.beginPath();
        for (const rect of blocks) { ctx.rect(rect.x, rect.y, rect.width, rect.height); }
        ctx.fillStrokeShape(shape);
      },
      fill: FEATURE_COLORS.block,
      strokeWidth: 0,
      listening: false
    }));
  }
  const byStatus = splitDotplotSnpEntriesByStatus(snps);
  drawDotplotSnpLines(layer, byStatus.rejected, FEATURE_COLORS.rejectedSnp, (ctx, snp) => {
    if (isX) { ctx.moveTo(snp.cx, snp.y0); ctx.lineTo(snp.cx, snp.y1); }
    else { ctx.moveTo(snp.x0, snp.cy); ctx.lineTo(snp.x1, snp.cy); }
  });
  drawDotplotSnpLines(layer, byStatus.pass, FEATURE_COLORS.snp, (ctx, snp) => {
    if (isX) { ctx.moveTo(snp.cx, snp.y0); ctx.lineTo(snp.cx, snp.y1); }
    else { ctx.moveTo(snp.x0, snp.cy); ctx.lineTo(snp.x1, snp.cy); }
  });
  layer.add(new Konva.Shape({
    sceneFunc(ctx, shape) {
      ctx.beginPath();
      drawRoundedRect(ctx, regionX, regionY, regionW, regionH, 2);
      ctx.fillStrokeShape(shape);
    },
    fillEnabled: false,
    stroke: "#000000",
    strokeWidth: 1,
    listening: false
  }));

  const hoverBlocks = blocks.map(rect => isX
    ? { x0: rect.x, x1: rect.x + rect.width, featureId: rect.featureId }
    : { y0: rect.y, y1: rect.y + rect.height, featureId: rect.featureId }
  );
  const hoverSnps = snps.map(entry => isX
    ? { cx: entry.cx, featureId: entry.featureId }
    : { cy: entry.cy, featureId: entry.featureId }
  );
  hoverBlocks.sort((a, b) => isX ? a.x0 - b.x0 : a.y0 - b.y0);
  hoverSnps.sort((a, b) => isX ? a.cx - b.cx : a.cy - b.cy);
  return { blocks: hoverBlocks, snps: hoverSnps, regionX, regionY, regionW, regionH, gffTracks: [] };
}

function drawDotplotGffTracks(layer, axis, localGeometry, sample, hoverIndex) {
  if (!sample) { return; }
  const tracks = getSampleGffTracks(sample);
  const isX = axis === "x";
  const queues = new Map();

  tracks.forEach((track, trackIndex) => {
    const color = getGffTrackColor(track.track_name);
    const trackOrigin = isX
      ? hoverIndex.regionY + hoverIndex.regionH + GFF_TRACK.topGap + trackIndex * (GFF_TRACK.height + GFF_TRACK.gap)
      : localGeometry.trackBoxX - localGeometry.yGffSideGap - trackIndex * (GFF_TRACK.height + GFF_TRACK.gap) - GFF_TRACK.height;
    const baseline = trackOrigin + GFF_TRACK.height / 2;
    const genes = [];

    layer.add(new Konva.Shape({
      sceneFunc(ctx) {
        ctx.save();
        ctx.strokeStyle = "#e5e7eb";
        ctx.lineWidth = 1;
        ctx.beginPath();
        if (isX) {
          ctx.moveTo(localGeometry.xZero, baseline);
          ctx.lineTo(localGeometry.xMax, baseline);
        } else {
          ctx.moveTo(baseline, localGeometry.yMaxPixel);
          ctx.lineTo(baseline, localGeometry.yZeroPixel);
        }
        ctx.stroke();
        ctx.restore();
      },
      listening: false
    }));

    for (const gene of track.features || []) {
      const p0 = isX
        ? mapXCoordinateToStagePx(gene.start_in_region, sample, localGeometry)
        : mapYCoordinateToStagePx(gene.start_in_region, sample, localGeometry);
      const p1 = isX
        ? mapXCoordinateToStagePx(gene.end_in_region, sample, localGeometry)
        : mapYCoordinateToStagePx(gene.end_in_region, sample, localGeometry);
      const along0 = Math.min(p0, p1);
      const length = Math.max(GFF_TRACK.minGeneWidthPx, Math.abs(p1 - p0));
      const rect = isX
        ? { x: along0, y: trackOrigin, width: length, height: GFF_TRACK.height }
        : { x: trackOrigin, y: along0, width: GFF_TRACK.height, height: length };
      if (!queues.has(color)) { queues.set(color, []); }
      queues.get(color).push(rect);
      genes.push({ x0: along0, x1eff: along0 + length, gene });
    }
    if (genes.length > 0) {
      genes.sort((a, b) => a.x0 - b.x0);
      hoverIndex.gffTracks.push(isX
        ? { y0: trackOrigin, y1: trackOrigin + GFF_TRACK.height, intervalIndex: buildGffGeneIntervalIndex(genes) }
        : { x0: trackOrigin, x1: trackOrigin + GFF_TRACK.height, intervalIndex: buildGffGeneIntervalIndex(genes) }
      );
    }
  });

  for (const [color, rects] of queues) {
    layer.add(new Konva.Shape({
      sceneFunc(ctx, shape) {
        ctx.globalAlpha = 0.85;
        ctx.beginPath();
        for (const rect of rects) { drawRoundedRect(ctx, rect.x, rect.y, rect.width, rect.height, 2); }
        ctx.globalAlpha = 1;
        ctx.fillStrokeShape(shape);
      },
      fill: color,
      strokeWidth: 0,
      listening: false
    }));
  }
}

function applyDotplotPaneLayout(geometry) {
  const container = document.querySelector(".dotplot-content");
  const surface = document.getElementById("dotplot-scroll-surface");
  const frame = document.querySelector(".dotplot-frame");
  const layout = document.getElementById("dotplot-ui-layout");
  if (!container || !surface || !frame || !layout) { return; }
  const values = {
    "--dotplot-surface-width": `${geometry.surfaceWidth}px`,
    "--dotplot-surface-height": `${geometry.surfaceHeight}px`,
    "--dotplot-y-pane-width": `${geometry.yContext.width}px`,
    "--dotplot-x-pane-height": `${geometry.xContext.height}px`,
    "--dotplot-matrix-width": `${geometry.matrixWidth}px`,
    "--dotplot-matrix-height": `${geometry.matrixHeight}px`,
    "--dotplot-matrix-viewport-height": `${geometry.matrixViewportHeight}px`,
    "--dotplot-viewport-width": `${geometry.scrollportWidth}px`,
    "--dotplot-top-breathing-space": `${DOTPLOT_FRAME_TOP_BREATHING_SPACE}px`,
    "--dotplot-gff-legend-height": `${GFF_LEGEND.height}px`
  };
  for (const [name, value] of Object.entries(values)) {
    layout.style.setProperty(name, value);
  }
  frame.style.height = `${geometry.scrollportHeight + DOTPLOT_FRAME_TOP_BREATHING_SPACE + 2}px`;
  surface.style.marginLeft = `${Math.max(0, (geometry.viewportWidth - geometry.scrollportWidth) / 2)}px`;
}

function getDotplotHorizontalScrollbarHeight(container) {
  // The scrollport has no border, so this difference is only the native
  // horizontal scrollbar. This runs on full redraws, never while scrolling.
  return Math.max(0, container.offsetHeight - container.clientHeight);
}

function renderDotplotGffLegend(xSampleData, ySampleData) {
  const legend = document.getElementById("dotplot-gff-legend");
  if (!legend) { return; }
  const hasSelectedGff = getSampleGffTracks(xSampleData).length > 0 || getSampleGffTracks(ySampleData).length > 0;
  const trackNames = getAllGffTrackNames();
  legend.innerHTML = "";
  legend.classList.toggle("hidden", !hasSelectedGff || trackNames.length === 0);
  if (!hasSelectedGff) { return; }
  for (const trackName of trackNames) {
    const item = document.createElement("span");
    item.className = "dotplot-gff-legend-item";
    const dot = document.createElement("span");
    dot.className = "dotplot-gff-legend-dot";
    dot.style.background = getGffTrackColor(trackName);
    item.append(dot, document.createTextNode(trackName));
    legend.appendChild(item);
  }
}

function redrawDotplotStage() {
  const img = document.getElementById("dotplot-svg-img");
  if (!img || !img.complete || img.naturalWidth === 0) {
    _currentDotplotGeometry = null;
    return;
  }
  initDotplotStage();
  const provisionalGeometry = computeDotplotGeometry();
  if (!provisionalGeometry) { _currentDotplotGeometry = null; return; }

  // First establish the frame dimensions, then measure the browser-owned
  // scrollbar once and recompute every dependent layout value coherently.
  applyDotplotPaneLayout(provisionalGeometry);
  const container = document.querySelector(".dotplot-content");
  const horizontalScrollbarHeight = container
    ? getDotplotHorizontalScrollbarHeight(container)
    : 0;
  const geometry = computeDotplotGeometry(horizontalScrollbarHeight);
  if (!geometry) { _currentDotplotGeometry = null; return; }
  _currentDotplotGeometry = geometry;

  const xSampleData = getSampleByName(_dotplotState.selectedX);
  const ySampleData = getSampleByName(_dotplotState.selectedY);
  applyDotplotPaneLayout(geometry);
  dotplotStage.size({ width: geometry.matrixWidth, height: geometry.matrixHeight });
  dotplotXStage.size({ width: geometry.xContext.width, height: geometry.xContext.height });
  dotplotYStage.size({ width: geometry.yContext.width, height: geometry.yContext.height });

  dotplotImageLayer.destroyChildren();
  dotplotImageLayer.add(new Konva.Image({
    x: geometry.imageX, y: geometry.imageY, image: img,
    width: geometry.imageWidth, height: geometry.imageHeight, listening: false
  }));

  dotplotXTrackLayer.destroyChildren();
  drawDotplotXAxis(dotplotXTrackLayer, geometry.xContext, xSampleData);
  const xHover = drawDotplotFeatureTrack(dotplotXTrackLayer, "x", geometry.xContext, xSampleData);
  drawDotplotGffTracks(dotplotXTrackLayer, "x", geometry.xContext, xSampleData, xHover);

  dotplotYTrackLayer.destroyChildren();
  drawDotplotYAxis(dotplotYTrackLayer, geometry.yContext, ySampleData);
  const yHover = drawDotplotFeatureTrack(dotplotYTrackLayer, "y", geometry.yContext, ySampleData);
  drawDotplotGffTracks(dotplotYTrackLayer, "y", geometry.yContext, ySampleData, yHover);
  _dotplotHoverIndex.xTrack = xHover;
  _dotplotHoverIndex.yTrack = yHover;
  _dotplotHoverIndexDirty = false;

  dotplotDebugLayer.destroyChildren();
  if (DOTPLOT_DEBUG_LAYOUT) {
    dotplotDebugLayer.add(new Konva.Shape({
      sceneFunc(ctx) {
        ctx.save();
        ctx.strokeStyle = DOTPLOT_TRACK.debugColor;
        ctx.lineWidth = DOTPLOT_TRACK.debugLineWidth;
        ctx.setLineDash([10, 4]);
        ctx.beginPath();
        ctx.moveTo(geometry.xZero, geometry.imageY); ctx.lineTo(geometry.xZero, geometry.imageY + geometry.imageHeight);
        ctx.moveTo(geometry.xMax, geometry.imageY); ctx.lineTo(geometry.xMax, geometry.imageY + geometry.imageHeight);
        ctx.moveTo(geometry.imageX, geometry.yZeroPixel); ctx.lineTo(geometry.imageX + geometry.imageWidth, geometry.yZeroPixel);
        ctx.moveTo(geometry.imageX, geometry.yMaxPixel); ctx.lineTo(geometry.imageX + geometry.imageWidth, geometry.yMaxPixel);
        ctx.stroke();
        ctx.restore();
      },
      listening: false
    }));
  }

  dotplotHighlightLayer.destroyChildren();
  dotplotHighlightLayer.add(_dotplotBlockIntersectionShape, _dotplotSnpProjectionShape);
  dotplotXHighlightLayer.destroyChildren();
  dotplotXHighlightLayer.add(_dotplotXBlockHighlightShape, _dotplotXSnpHighlightShape);
  dotplotYHighlightLayer.destroyChildren();
  dotplotYHighlightLayer.add(_dotplotYBlockHighlightShape, _dotplotYSnpHighlightShape);
  updateDotplotHighlightShapes();
  renderDotplotGffLegend(xSampleData, ySampleData);

  for (const viewerId of ["dotplot-viewer", "dotplot-x-viewer", "dotplot-y-viewer"]) {
    document.getElementById(viewerId)?.classList.remove("hidden");
  }
  dotplotStage.draw();
  dotplotXStage.draw();
  dotplotYStage.draw();
  restorePendingDotplotViewportCenter();
}

// Binary search: returns index of first element where arr[i].cx >= target.
function lowerBoundDotplotCx(arr, target) {
  let lo = 0;
  let hi = arr.length;
  while (lo < hi) {
    const mid = (lo + hi) >>> 1;
    if (arr[mid].cx < target) { lo = mid + 1; } else { hi = mid; }
  }
  return lo;
}

// Binary search: returns index of first element where arr[i].cy >= target.
function lowerBoundDotplotCy(arr, target) {
  let lo = 0;
  let hi = arr.length;
  while (lo < hi) {
    const mid = (lo + hi) >>> 1;
    if (arr[mid].cy < target) { lo = mid + 1; } else { hi = mid; }
  }
  return lo;
}

// Binary search: returns index of first element where arr[i].x0 >= target.
function lowerBoundDotplotX0(arr, target) {
  let lo = 0;
  let hi = arr.length;
  while (lo < hi) {
    const mid = (lo + hi) >>> 1;
    if (arr[mid].x0 < target) { lo = mid + 1; } else { hi = mid; }
  }
  return lo;
}

// Binary search: returns index of first element where arr[i].y0 >= target.
function lowerBoundDotplotY0(arr, target) {
  let lo = 0;
  let hi = arr.length;
  while (lo < hi) {
    const mid = (lo + hi) >>> 1;
    if (arr[mid].y0 < target) { lo = mid + 1; } else { hi = mid; }
  }
  return lo;
}

// Resolves which feature (if any) the pointer is hovering over in one local
// context stage.  The caller supplies the stage axis, so no sticky-offset
// compensation is needed.
// SNPs are prioritized over blocks (closer distance wins; tolerance = SNP_POINTER_TOLERANCE_PX).
// Returns { featureType: "block"|"snp", featureId } or null.
function resolveDotplotHoveredFeature(axis, pointerX, pointerY) {
  if (_dotplotHoverIndexDirty) {
    return null;
  }

  const track = axis === "x" ? _dotplotHoverIndex.xTrack : _dotplotHoverIndex.yTrack;

  if (
    track.regionH > 0 &&
    pointerX >= track.regionX && pointerX <= track.regionX + track.regionW &&
    pointerY >= track.regionY && pointerY <= track.regionY + track.regionH
  ) {
    if (axis === "x") {
    // SNPs: tolerance scan around pointerX.
    const lo = lowerBoundDotplotCx(track.snps, pointerX - SNP_POINTER_TOLERANCE_PX);
    let closestSnpId = null;
    let closestDist  = SNP_POINTER_TOLERANCE_PX + 1;
    for (let j = lo; j < track.snps.length; j++) {
      const s = track.snps[j];
      if (s.cx > pointerX + SNP_POINTER_TOLERANCE_PX) { break; }
      const d = Math.abs(pointerX - s.cx);
      if (d < closestDist) { closestDist = d; closestSnpId = s.featureId; }
    }
    if (closestSnpId !== null) {
      return { featureType: "snp", featureId: closestSnpId };
    }
    // Blocks: last block whose x0 <= pointerX, check x1.
    const bi = lowerBoundDotplotX0(track.blocks, pointerX) - 1;
    if (bi >= 0 && pointerX <= track.blocks[bi].x1) {
      return { featureType: "block", featureId: track.blocks[bi].featureId };
    }
    return null;
    }
    // SNPs: tolerance scan around pointerY.
    const lo = lowerBoundDotplotCy(track.snps, pointerY - SNP_POINTER_TOLERANCE_PX);
    let closestSnpId = null;
    let closestDist  = SNP_POINTER_TOLERANCE_PX + 1;
    for (let j = lo; j < track.snps.length; j++) {
      const s = track.snps[j];
      if (s.cy > pointerY + SNP_POINTER_TOLERANCE_PX) { break; }
      const d = Math.abs(pointerY - s.cy);
      if (d < closestDist) { closestDist = d; closestSnpId = s.featureId; }
    }
    if (closestSnpId !== null) {
      return { featureType: "snp", featureId: closestSnpId };
    }
    // Blocks: last block whose y0 <= pointerY, check y1.
    const bi = lowerBoundDotplotY0(track.blocks, pointerY) - 1;
    if (bi >= 0 && pointerY <= track.blocks[bi].y1) {
      return { featureType: "block", featureId: track.blocks[bi].featureId };
    }
    return null;
  }

  return null;
}

function resolveDotplotGffGene(axis, pointerX, pointerY) {
  if (_dotplotHoverIndexDirty) { return null; }
  const tracks = (axis === "x" ? _dotplotHoverIndex.xTrack : _dotplotHoverIndex.yTrack).gffTracks || [];
  for (const track of tracks) {
    const isWithinTrack = axis === "x"
      ? pointerY >= track.y0 && pointerY <= track.y1
      : pointerX >= track.x0 && pointerX <= track.x1;
    if (isWithinTrack) {
      return findGffGeneAtX(track.intervalIndex, axis === "x" ? pointerX : pointerY);
    }
  }
  return null;
}

function showDotplotGffGeneTooltip(gene, stage, pointerX, pointerY) {
  const tooltip = document.getElementById("gff-gene-tooltip");
  if (!tooltip) { return; }
  tooltip.textContent = formatGffGeneAttributes(gene.attributes);
  tooltip.style.display = "block";
  const stageRect = stage.container().getBoundingClientRect();
  const tooltipRect = tooltip.getBoundingClientRect();
  const margin = 8;
  const left = clampValue(stageRect.left + pointerX + 12, margin, window.innerWidth - tooltipRect.width - margin);
  const top = clampValue(stageRect.top + pointerY + 12, margin, window.innerHeight - tooltipRect.height - margin);
  tooltip.style.left = `${left}px`;
  tooltip.style.top = `${top}px`;
}

// Applies a resolved dotplot hover, guarded by a key to avoid redundant sidebar updates.
function applyDotplotResolvedHover(resolved) {
  const key = resolved ? `${resolved.featureType}:${resolved.featureId}` : null;
  if (key === _lastResolvedDotplotHoverKey) {
    return;
  }
  _lastResolvedDotplotHoverKey = key;

  if (!resolved) {
    clearHoveredFeature();
    return;
  }

  setHoveredFeature(resolved.featureType, resolved.featureId);
}

// Computes independent local-stage highlight geometry for both genomic panes.
function getDotplotHighlightGeometries(featureType, featureId) {
  const x = { blockGeoms: [], snpGeoms: [] };
  const y = { blockGeoms: [], snpGeoms: [] };

  if (!featureId || _dotplotHoverIndexDirty) {
    return { x, y };
  }

  const { xTrack, yTrack } = _dotplotHoverIndex;

  if (featureType === "block") {
    // Block highlights span the full track cross-axis with TRACK_HIGHLIGHT_INSET,
    // matching browser-mode getBlockHighlightGeometries (0.5 px inset from region border).
    for (const b of xTrack.blocks) {
      if (b.featureId === featureId) {
        x.blockGeoms.push({
          x:      b.x0,
          y:      xTrack.regionY + TRACK_HIGHLIGHT_INSET,
          width:  b.x1 - b.x0,
          height: xTrack.regionH - 2 * TRACK_HIGHLIGHT_INSET
        });
      }
    }
    for (const b of yTrack.blocks) {
      if (b.featureId === featureId) {
        y.blockGeoms.push({
          x:      yTrack.regionX + TRACK_HIGHLIGHT_INSET,
          y:      b.y0,
          width:  yTrack.regionW - 2 * TRACK_HIGHLIGHT_INSET,
          height: b.y1 - b.y0
        });
      }
    }
  } else if (featureType === "snp") {
    // SNP highlights span the same inset as feature lines.
    for (const s of xTrack.snps) {
      if (s.featureId === featureId) {
        x.snpGeoms.push({
          cx: s.cx,
          y0: xTrack.regionY + TRACK_GEOMETRY.featureInset,
          y1: xTrack.regionY + xTrack.regionH - TRACK_GEOMETRY.featureInset
        });
      }
    }
    for (const s of yTrack.snps) {
      if (s.featureId === featureId) {
        y.snpGeoms.push({
          cy: s.cy,
          x0: yTrack.regionX + TRACK_GEOMETRY.featureInset,
          x1: yTrack.regionX + yTrack.regionW - TRACK_GEOMETRY.featureInset
        });
      }
    }
  }

  return { x, y };
}

// Generic value clamp helper.
function clampValue(value, min, max) {
  return Math.max(min, Math.min(max, value));
}

// Resolves the block feature ID to use for the dotplot intersection overlay.
// For a pinned block, returns its featureId directly.
// For a pinned SNP, looks up block_id from featureGroups then resolves it to
// a block featureId via searchIndexes.blockIdToFeatureId.
function _resolveDotplotCenterBlockFeatureId() {
  if (!state.pinnedFeatureId || !state.pinnedFeatureType) { return null; }
  if (state.pinnedFeatureType === "block") {
    return state.pinnedFeatureId;
  }
  if (state.pinnedFeatureType === "snp") {
    const entries = state.featureGroups.get(state.pinnedFeatureId);
    if (!entries || entries.length === 0) { return null; }
    const blockId = entries[0].info.block_id;
    if (blockId == null) { return null; }
    return searchIndexes.blockIdToFeatureId.get(Number(blockId)) ?? null;
  }
  return null;
}


function resetActiveViewerZoom() {
  if (isDotplotModeActive()) {
    _dotplotViewportRestoreEpoch += 1;
    _pendingDotplotViewportRestore = null;
    _dotplotState.zoom = 1;

    redrawDotplotStage();

    requestAnimationFrame(() => {
      const container = document.querySelector(".dotplot-content");
      if (container) {
        container.scrollLeft = 0;
        container.scrollTop = 0;
      }
    });

    return;
  }

  state.zoomX = getInitialZoomX();
  state.scrollX = 0;
  redrawStage();
}


function moveDotplotByViewportFraction(direction, fraction = 0.1) {
  const container = document.querySelector(".dotplot-content");
  if (!container) {
    return;
  }

  const stepPx = (_currentDotplotGeometry?.matrixViewportWidth || container.clientWidth) * fraction;
  container.scrollLeft = clampValue(
    container.scrollLeft + direction * stepPx,
    0,
    container.scrollWidth - container.clientWidth
  );
}

// Centers the dotplot viewport on the pinned feature's block intersection
// rectangle, adapting the zoom level so the rectangle occupies a meaningful
// fraction of the visible viewport.
//
// Algorithm:
//   1. Compute the intersection rect at current zoom to get its natural size.
//   2. Pick the zoom that makes rect cover TARGET_COVERAGE of the viewport,
//      independently for width and height; take the smaller (zoom-in less).
//   3. Clamp to [DOTPLOT_ZOOM_MIN, DOTPLOT_ZOOM_MAX].
//   4. Apply zoom, redraw synchronously, recompute rect.
//   5. Scroll both axes so the rect center is in the middle of the viewport.
function centerDotplotOnPinnedFeature() {
  if (!isDotplotModeActive()) { return; }

  const blockFeatureId = _resolveDotplotCenterBlockFeatureId();
  if (!blockFeatureId) { return; }

  // Step 1 — compute rect at current zoom to know its real-unit size.
  const rectAtCurrentZoom = _computeDotplotBlockIntersection(blockFeatureId);
  if (!rectAtCurrentZoom) { return; }

  const container = document.querySelector(".dotplot-content");
  if (!container) { return; }

  // Step 2 — adaptive zoom: use only the usable matrix viewport; frozen
  // context panes are excluded from the target coverage calculation.
  const TARGET_COVERAGE = 0.50;
  const viewW = _currentDotplotGeometry?.matrixViewportWidth || container.clientWidth;
  const viewH = _currentDotplotGeometry?.matrixViewportHeight || container.clientHeight;
  const currentRectW = rectAtCurrentZoom.width;
  const currentRectH = rectAtCurrentZoom.height;

  // Desired zoom independently for each axis, then take the smaller one so
  // the entire rect fits while still being as large as possible.
  const desiredZoomX = (currentRectW > 0 && viewW > 0)
    ? (_dotplotState.zoom * viewW * TARGET_COVERAGE) / currentRectW
    : DOTPLOT_ZOOM_MAX;
  const desiredZoomY = (currentRectH > 0 && viewH > 0)
    ? (_dotplotState.zoom * viewH * TARGET_COVERAGE) / currentRectH
    : DOTPLOT_ZOOM_MAX;

  let targetZoom = Math.min(desiredZoomX, desiredZoomY);
  targetZoom = clampValue(targetZoom, DOTPLOT_ZOOM_MIN, DOTPLOT_ZOOM_MAX);

  // Step 3 — apply zoom and redraw synchronously so geometry is up to date.
  _dotplotState.zoom = targetZoom;
  redrawDotplotStage();

  // Step 4 — recompute rect with the new zoom.
  const rect = _computeDotplotBlockIntersection(blockFeatureId);
  if (!rect) { return; }

  // Step 5 — apply scroll in the next frame so the browser has committed the
  // new canvas dimensions.
  requestAnimationFrame(() => {
    const rectCenterX = rect.x + rect.width / 2;
    const rectCenterY = rect.y + rect.height / 2;

    const matrixViewportWidth = _currentDotplotGeometry?.matrixViewportWidth || container.clientWidth;
    const matrixViewportHeight = _currentDotplotGeometry?.matrixViewportHeight || container.clientHeight;
    container.scrollLeft = clampValue(
      rectCenterX - matrixViewportWidth / 2,
      0,
      container.scrollWidth - container.clientWidth
    );

    container.scrollTop = clampValue(
      rectCenterY - matrixViewportHeight / 2,
      0,
      container.scrollHeight - container.clientHeight
    );

    const containerRect = container.getBoundingClientRect();
    const containerCenterYInPage = window.scrollY + containerRect.top + containerRect.height / 2;
    const targetWindowScrollY = containerCenterYInPage - window.innerHeight / 2;

    window.scrollTo({
      top: clampValue(
        targetWindowScrollY,
        0,
        document.documentElement.scrollHeight - window.innerHeight
      ),
      behavior: "smooth"
    });

    showToast("R: reset zoom · F: center · Tab/Shift+Tab: navigate");
  });
}

function centerPinnedFeature() {
  if (!state.pinnedFeatureId || !state.pinnedFeatureType) {
    return;
  }

  if (isDotplotModeActive()) {
    centerDotplotOnPinnedFeature();
    return;
  }

  const range = searchIndexes.featureIdToRegionRange.get(state.pinnedFeatureId);
  if (!range) {
    return;
  }

  centerRegionOnRange(range.start, range.end);
}

function _getDotplotMappedBlockBounds(featureId) {
  const geometry = _currentDotplotGeometry;
  if (!geometry) { return null; }

  const xSampleData = getSampleByName(_dotplotState.selectedX);
  const ySampleData = getSampleByName(_dotplotState.selectedY);
  if (!xSampleData || !ySampleData) { return null; }

  const xBlock = xSampleData.blocks.find(b => b.feature_id === featureId);
  const yBlock = ySampleData.blocks.find(b => b.feature_id === featureId);
  if (!xBlock || !yBlock) { return null; }

  const xLeft = mapXCoordinateToStagePx(xBlock.block_start_in_region, xSampleData, geometry);
  const xRight = mapXCoordinateToStagePx(xBlock.block_end_in_region, xSampleData, geometry);

  const yStart = mapYCoordinateToStagePx(yBlock.block_start_in_region, ySampleData, geometry);
  const yEnd = mapYCoordinateToStagePx(yBlock.block_end_in_region, ySampleData, geometry);

  const x0 = Math.min(xLeft, xRight);
  const x1 = Math.max(xLeft, xRight);
  const y0 = Math.min(yStart, yEnd);
  const y1 = Math.max(yStart, yEnd);

  return { geometry, x0, x1, y0, y1 };
}

function _computeDotplotBlockIntersection(featureId) {
  const bounds = _getDotplotMappedBlockBounds(featureId);
  if (!bounds) { return null; }

  const { x0, x1, y0, y1 } = bounds;

  return {
    x: x0,
    y: y0,
    width: Math.max(DOTPLOT_RENDERING.intersectionMinSizePx, x1 - x0),
    height: Math.max(DOTPLOT_RENDERING.intersectionMinSizePx, y1 - y0)
  };
}

// Computes the stage-space rectangle for the block projection overlay:
// the area on the SVG image that corresponds to the given block's coordinate
// interval on both the selected X and Y samples.
// Returns { x, y, width, height } in stage pixels, or null.
function _computeDotplotBlockProjection(featureId) {
  const bounds = _getDotplotMappedBlockBounds(featureId);
  if (!bounds) { return null; }

  const { geometry, x0, x1, y0, y1 } = bounds;

  return {
    vertical: {
      x: x0,
      y: y0,
      width: Math.max(DOTPLOT_RENDERING.intersectionMinSizePx, x1 - x0),
      height: geometry.yZeroPixel - y0
    },
    horizontal: {
      x: geometry.xZero,
      y: y0,
      width: x1 - geometry.xZero,
      height: Math.max(DOTPLOT_RENDERING.intersectionMinSizePx, y1 - y0)
    }
  };
}


function _computeDotplotSnpProjection(featureId) {
  if (!shouldDisplaySnp(featureId)) {
    return null;
  }

  const geometry = _currentDotplotGeometry;
  if (!geometry) { return null; }

  const xSampleData = getSampleByName(_dotplotState.selectedX);
  const ySampleData = getSampleByName(_dotplotState.selectedY);
  if (!xSampleData || !ySampleData) { return null; }

  const xSnp = xSampleData.snps.find(s => s.feature_id === featureId);
  const ySnp = ySampleData.snps.find(s => s.feature_id === featureId);
  if (!xSnp || !ySnp) { return null; }

  const x = mapXCoordinateToStagePx(xSnp.pos_in_region, xSampleData, geometry);
  const y = mapYCoordinateToStagePx(ySnp.pos_in_region, ySampleData, geometry);

  return {
    x,
    y,
    xZero: geometry.xZero,
    yZeroPixel: geometry.yZeroPixel
  };
}

// Updates matrix projections and both local context highlight layers.  This is
// triggered by feature interaction or a full redraw, never by native scrolling.
function updateDotplotHighlightShapes() {
  if (
    !dotplotHighlightLayer || !dotplotXHighlightLayer || !dotplotYHighlightLayer ||
    !_dotplotXBlockHighlightShape || !_dotplotXSnpHighlightShape ||
    !_dotplotYBlockHighlightShape || !_dotplotYSnpHighlightShape
  ) {
    return;
  }

  const displayed = getDisplayedFeature();
  const color = displayed && displayed.source === "pin"
    ? FEATURE_COLORS.highlightPinned
    : FEATURE_COLORS.highlightHover;

  let xGeoms = { blockGeoms: [], snpGeoms: [] };
  let yGeoms = { blockGeoms: [], snpGeoms: [] };

  if (
    displayed &&
    !(displayed.featureType === "snp" && !shouldDisplaySnp(displayed.featureId))
  ) {
    const result = getDotplotHighlightGeometries(displayed.featureType, displayed.featureId);
    xGeoms = result.x;
    yGeoms = result.y;
  }

  _dotplotXBlockHighlightColor = color;
  _dotplotXBlockHighlightGeoms = xGeoms.blockGeoms;
  _dotplotXBlockHighlightShape.visible(xGeoms.blockGeoms.length > 0);
  _dotplotXSnpHighlightColor = color;
  _dotplotXSnpHighlightGeoms = xGeoms.snpGeoms;
  _dotplotXSnpHighlightShape.visible(xGeoms.snpGeoms.length > 0);
  _dotplotYBlockHighlightColor = color;
  _dotplotYBlockHighlightGeoms = yGeoms.blockGeoms;
  _dotplotYBlockHighlightShape.visible(yGeoms.blockGeoms.length > 0);
  _dotplotYSnpHighlightColor = color;
  _dotplotYSnpHighlightGeoms = yGeoms.snpGeoms;
  _dotplotYSnpHighlightShape.visible(yGeoms.snpGeoms.length > 0);

  // Block intersection overlay: shown only for block features in dotplot mode.
  if (_dotplotBlockIntersectionShape) {
    if (displayed && displayed.featureType === "block") {
      _dotplotBlockIntersectionGeom = _computeDotplotBlockProjection(displayed.featureId);
    } else {
      _dotplotBlockIntersectionGeom = null;
    }
    _dotplotBlockIntersectionShape.visible(_dotplotBlockIntersectionGeom !== null);
  }

  if (_dotplotSnpProjectionShape) {
    if (
      displayed &&
      displayed.featureType === "snp" &&
      shouldDisplaySnp(displayed.featureId)
    ) {
      _dotplotSnpProjectionGeom = _computeDotplotSnpProjection(displayed.featureId);
    } else {
      _dotplotSnpProjectionGeom = null;
    }
    _dotplotSnpProjectionShape.visible(_dotplotSnpProjectionGeom !== null);
  }

  dotplotHighlightLayer.batchDraw();
  dotplotXHighlightLayer.batchDraw();
  dotplotYHighlightLayer.batchDraw();
}

function renderDotplot() {
  const img = document.getElementById("dotplot-svg-img");

  if (!img) {
    return;
  }

  _currentDotplotGeometry = null;

  const pair = findDotplotPair(_dotplotState.selectedY, _dotplotState.selectedX);

  if (!pair) {
    updateDotplotStatusMessage("Select a compatible X/Y sample.");
    return;
  }

  _dotplotHoverIndexDirty = true;
  _lastResolvedDotplotHoverKey = null;

  const pairKey = getSelectedDotplotPairKey();
  _renderedDotplotPairKey = pairKey;
  img.onload = () => {
    if (_renderedDotplotPairKey === pairKey) {
      requestDotplotRedraw();
    }
  };
  img.src = pair.svg_rel_path;

  updateDotplotStatusMessage("");

  if (img.complete && img.naturalWidth > 0) {
    requestDotplotRedraw();
  }
}

function isCompleteValidDotplotPair() {
  return Boolean(findDotplotPair(_dotplotState.selectedY, _dotplotState.selectedX));
}

function getCompatibleDotplotXSamples(ySample) {
  return dotplotPairIndexes.validPairsByY.get(ySample) || new Set();
}

function getCompatibleDotplotYSamples(xSample) {
  return dotplotPairIndexes.validPairsByX.get(xSample) || new Set();
}

function getDotplotLabelState(axis, sampleName) {
  const isX = axis === "x";
  const selectedValue = isX ? _dotplotState.selectedX : _dotplotState.selectedY;
  const oppositeValue = isX ? _dotplotState.selectedY : _dotplotState.selectedX;

  const isSelected = selectedValue === sampleName;
  const isPending = _dotplotState.pendingAxis === axis && isSelected;

  if (isPending) {
    return { enabled: true, className: "pending" };
  }

  if (isCompleteValidDotplotPair()) {
    return {
      enabled: true,
      className: isSelected ? `active-${axis}` : ""
    };
  }

  if (_dotplotState.pendingAxis === "x" && !isX) {
    const compatibleY = getCompatibleDotplotYSamples(_dotplotState.selectedX);
    const enabled = compatibleY.has(sampleName);
    return { enabled, className: enabled ? "" : "disabled" };
  }

  if (_dotplotState.pendingAxis === "y" && isX) {
    const compatibleX = getCompatibleDotplotXSamples(_dotplotState.selectedY);
    const enabled = compatibleX.has(sampleName);
    return { enabled, className: enabled ? "" : "disabled" };
  }

  return {
    enabled: true,
    className: isSelected ? `active-${axis}` : ""
  };
}

function createDotplotSampleLabel(axis, sampleName) {
  const labelState = getDotplotLabelState(axis, sampleName);
  const button = document.createElement("button");

  button.type = "button";
  button.textContent = sampleName;
  button.dataset.axis = axis;
  button.dataset.sample = sampleName;
  button.className = `dotplot-sample-label ${labelState.className}`.trim();
  button.disabled = !labelState.enabled;

  button.addEventListener("click", () => {
    handleDotplotSampleLabelClick(axis, sampleName);
  });

  return button;
}

function renderDotplotSampleLabels() {
  const xContainer = document.getElementById("dotplot-x-labels");
  const yContainer = document.getElementById("dotplot-y-labels");

  if (!xContainer || !yContainer) {
    return;
  }

  xContainer.innerHTML = "";
  yContainer.innerHTML = "";

  for (const sampleName of dotplotPairIndexes.samples) {
    xContainer.appendChild(createDotplotSampleLabel("x", sampleName));
    yContainer.appendChild(createDotplotSampleLabel("y", sampleName));
  }
}

function handleDotplotSampleLabelClick(axis, sampleName) {
  if (axis === "x") {
    handleDotplotXLabelClick(sampleName);
  } else {
    handleDotplotYLabelClick(sampleName);
  }

  renderDotplotSampleLabels();
  updateDotplotPendingPairVisualState();
}

function handleDotplotXLabelClick(sampleName) {
  if (isCompleteValidDotplotPair() && _dotplotState.selectedX === sampleName) {
    return;
  }

  if (_dotplotState.selectedY === null) {
    _dotplotState.selectedX = sampleName;
    _dotplotState.pendingAxis = "x";
    updateDotplotStatusMessage("Select a compatible Y sample.");
    return;
  }

  if (_dotplotState.pendingAxis === "y") {
    const pair = findDotplotPair(_dotplotState.selectedY, sampleName);
    if (!pair) {
      return;
    }

    _dotplotState.selectedX = sampleName;
    _dotplotState.pendingAxis = null;
    loadSelectedDotplotPair();
    return;
  }

  _dotplotState.selectedX = sampleName;
  _dotplotState.selectedY = null;
  _dotplotState.pendingAxis = "x";
  updateDotplotStatusMessage("Select a compatible Y sample.");
}

function handleDotplotYLabelClick(sampleName) {
  if (isCompleteValidDotplotPair() && _dotplotState.selectedY === sampleName) {
    return;
  }

  if (_dotplotState.selectedX === null) {
    _dotplotState.selectedY = sampleName;
    _dotplotState.pendingAxis = "y";
    updateDotplotStatusMessage("Select a compatible X sample.");
    return;
  }

  if (_dotplotState.pendingAxis === "x") {
    const pair = findDotplotPair(sampleName, _dotplotState.selectedX);
    if (!pair) {
      return;
    }

    _dotplotState.selectedY = sampleName;
    _dotplotState.pendingAxis = null;
    loadSelectedDotplotPair();
    return;
  }

  _dotplotState.selectedY = sampleName;
  _dotplotState.selectedX = null;
  _dotplotState.pendingAxis = "y";
  updateDotplotStatusMessage("Select a compatible X sample.");
}

function updateDotplotStatusMessage(message) {
  const msg = document.getElementById("dotplot-status-msg");
  if (!msg) {
    return;
  }

  msg.textContent = message;
  msg.classList.toggle("hidden", message === "");
}

function updateDotplotPendingPairVisualState() {
  const container = document.querySelector(".dotplot-content");
  if (!container) {
    return;
  }

  container.classList.toggle(
    "is-pending-pair",
    !isCompleteValidDotplotPair()
  );
}

function loadSelectedDotplotPair() {
  const pairKey = getSelectedDotplotPairKey();
  const previousViewportCenter = pairKey && pairKey !== _renderedDotplotPairKey
    ? captureDotplotViewportCenter()
    : null;
  const epoch = ++_dotplotViewportRestoreEpoch;
  _pendingDotplotViewportRestore = previousViewportCenter
    ? { ...previousViewportCenter, pairKey, epoch }
    : null;
  updateDotplotPendingPairVisualState();
  renderDotplot();
}

function setupDotplotUI() {
  const pairs = getDotplotPairs();
  const controls = document.querySelector(".dotplot-controls");

  if (pairs.length === 0) {
    if (controls) {
      controls.classList.add("hidden");
    }
    updateDotplotStatusMessage("No dotplots are available.");
    return;
  }

  const firstPair = pairs[0];
  _dotplotState.selectedY = firstPair.y_sample;
  _dotplotState.selectedX = firstPair.x_sample;
  _dotplotState.pendingAxis = null;

  renderDotplotSampleLabels();
  renderDotplot();
}

function setViewerMode(mode) {
  const viewerCanvas       = document.getElementById("viewer");
  const viewerToolbar      = document.querySelector(".viewer-toolbar");
  const alignmentPanel     = document.getElementById("alignment-panel");
  const dotplotPanel       = document.getElementById("dotplot-panel");
  const browserBtn         = document.getElementById("browser-mode-btn");
  const dotplotBtn         = document.getElementById("dotplot-mode-btn");

  const isBrowser = mode === "browser";
  
  updateViewerModeInfoTooltip(mode);

  if (viewerCanvas) {
    viewerCanvas.classList.toggle("hidden", !isBrowser);
  }
  if (dotplotPanel) {
    dotplotPanel.classList.toggle("hidden", isBrowser);
  }
  if (browserBtn) {
    browserBtn.classList.toggle("active", isBrowser);
  }
  if (dotplotBtn) {
    dotplotBtn.classList.toggle("active", !isBrowser);
  }

  updateSearchModeAvailability(isBrowser);

  if (isBrowser) {
    requestStageRedraw();
    requestActiveAlignmentViewerUpdate();
  } else {
    // Give the panel time to become visible before computing image dimensions.
    requestDotplotRedraw();
    requestActiveAlignmentViewerUpdate();
  }
  updateFeatureNavigationButtons();
  syncSidebarHeightToViewerColumn();
}

function setupModeSwitch() {
  const browserBtn = document.getElementById("browser-mode-btn");
  const dotplotBtn = document.getElementById("dotplot-mode-btn");

  if (browserBtn) {
    browserBtn.addEventListener("click", () => {
      setViewerMode("browser");
    });
  }
  if (dotplotBtn) {
    dotplotBtn.addEventListener("click", () => {
      setViewerMode("dotplot");
    });
  }
}
