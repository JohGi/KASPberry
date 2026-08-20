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

function getDotplotYSamples() {
  const seen = new Set();
  const result = [];
  for (const pair of getDotplotPairs()) {
    if (!seen.has(pair.y_sample)) {
      seen.add(pair.y_sample);
      result.push(pair.y_sample);
    }
  }
  return result;
}

function getDotplotXSamplesForY(ySample) {
  return getDotplotPairs()
    .filter(pair => pair.y_sample === ySample)
    .map(pair => pair.x_sample);
}

function findDotplotPair(ySample, xSample) {
  if (!ySample || !xSample) {
    return null;
  }

  return dotplotPairIndexes.pairByYX.get(`${ySample}::${xSample}`) || null;
}

function getSelectedDotplotPair() {
  return findDotplotPair(_dotplotState.selectedY, _dotplotState.selectedX);
}

function getSampleByName(sampleName) {
  return REGION_DATA.samples.find(s => s.sample === sampleName) || null;
}

// Returns the displayed image size, computed from naturalWidth/Height.
// Applies zoom, then clamps to available container space to avoid overflow at zoom=1.
// At zoom > 1 the image may exceed the container and the panel will scroll.
function getDotplotImageDisplaySize() {
  const img = document.getElementById("dotplot-svg-img");
  if (!img || !img.complete || img.naturalWidth === 0) {
    return null;
  }
  const container = document.querySelector(".dotplot-content");
  const containerPadding = 24; // 12 px each side
  // Use full available container width minus the y-track gutter and track gap.
  const availContainerW = container ? container.clientWidth - containerPadding : 800;
  
  const reservedAxisSpace = 80;

  const maxW = Math.max(
    100,
    availContainerW - DOTPLOT_TRACK.yTrackWidth - reservedAxisSpace
  );

  const maxH = Math.max(
    100,
    Math.floor(window.innerHeight * 0.9) - DOTPLOT_TRACK.xTrackHeight - reservedAxisSpace
  );

  // Base size: image scaled to fit inside maxW × maxH while preserving aspect ratio.
  let w = img.naturalWidth;
  let h = img.naturalHeight;

  if (w > maxW) {
    h = Math.round(h * maxW / w);
    w = maxW;
  }
  if (h > maxH) {
    w = Math.round(w * maxH / h);
    h = maxH;
  }

  // Apply zoom on top of the fitted base size.
  w = Math.round(w * _dotplotState.zoom);
  h = Math.round(h * _dotplotState.zoom);

  return { imageWidth: Math.max(1, w), imageHeight: Math.max(1, h) };
}

// Returns the total pixel width that Y-sample GFF tracks occupy to the left of the Y region,
// including the topGap used as a side-gap between GFF tracks and the Y region border.
function getDotplotYGffTotalWidth(ySampleData) {
  // Returns the total width for Y-sample GFF tracks, not including the side gap.
  if (!ySampleData) { return 0; }
  const n = getSampleGffTracks(ySampleData).length;
  if (n === 0) { return 0; }
  return n * (GFF_TRACK.height + GFF_TRACK.gap);
}

// Returns the total pixel height that X-sample GFF tracks occupy below the X region,
// including topGap and an optional legend row.
function getDotplotXGffTotalHeight(xSampleData) {
  const trackCount = xSampleData ? getSampleGffTracks(xSampleData).length : 0;
  const legendTopGap = getAllGffTrackNames().length > 0 ? 20 : 0;
  const legendH = getAllGffTrackNames().length > 0
    ? legendTopGap + GFF_LEGEND.height
    : 0;
  if (trackCount === 0 && legendH === 0) { return 0; }
  return (trackCount > 0 ? GFF_TRACK.topGap + trackCount * (GFF_TRACK.height + GFF_TRACK.gap) : 0)
    + legendH;
}

// Computes the full geometry for the dotplot Konva stage.
// All coordinates are in stage space:
//   y-track region:  x = [TRACK_FEATURE_INSET, yTrackWidth-TRACK_FEATURE_INSET],
//                  y = [yMaxPixel, yZeroPixel].
//   image occupies x = [yTrackWidth+DOTPLOT_TRACK_GAP, …],  y = [0, imageHeight].
//   x-track region:  x = [xZero, xMax],
//                  y = [imageHeight+DOTPLOT_TRACK_GAP+TRACK_FEATURE_INSET, …].
// The axis-bounds ratios (DOTPLOT_AXIS_BOUNDS) are applied to imageWidth/Height so
// coordinate mapping is always relative to the image, regardless of gap size.
function computeDotplotGeometry() {
  const size = getDotplotImageDisplaySize();
  if (!size) {
    return null;
  }
  const { imageWidth, imageHeight } = size;
  const { yTrackWidth, xTrackHeight } = DOTPLOT_TRACK;

  const xSampleData = getSampleByName(_dotplotState.selectedX);
  const ySampleData = getSampleByName(_dotplotState.selectedY);

  // Extra horizontal space on the left for Y-sample GFF tracks (not including side gap).
  const yGffWidth = getDotplotYGffTotalWidth(ySampleData);
  // Side gap between Y region and GFF tracks (same as GFF_TRACK.topGap for symmetry).
  const yGffSideGap = yGffWidth > 0 ? GFF_TRACK.topGap : 0;
  // Extra vertical space below the X region for X-sample GFF tracks + legend.
  const xGffHeight = getDotplotXGffTotalHeight(xSampleData);

  const xAxisGap = getDotplotXAxisGap();
  const yAxisGap = getDotplotYAxisGap(ySampleData, imageHeight);

  // Image is offset right by the y-track width + y-GFF gutter + side gap + gap.
  const imageX = yGffWidth + yGffSideGap + yTrackWidth + yAxisGap;
  const imageY = 0;

  const xZero     = imageX + imageWidth  * DOTPLOT_AXIS_BOUNDS.xZeroRatio;
  const xMax      = imageX + imageWidth  * DOTPLOT_AXIS_BOUNDS.xMaxRatio;
  const yZeroPixel = imageY + imageHeight * (1 - DOTPLOT_AXIS_BOUNDS.yZeroRatio);
  const yMaxPixel  = imageY + imageHeight * (1 - DOTPLOT_AXIS_BOUNDS.yMaxRatio);

  const stageStrokePadding = 1;

  return {
    stageWidth:  yGffWidth + yGffSideGap + yTrackWidth + yAxisGap + imageWidth + stageStrokePadding,
    stageHeight: imageHeight + xAxisGap + xTrackHeight + xGffHeight + stageStrokePadding,
    imageX,
    imageY,
    imageWidth,
    imageHeight,
    yTrackWidth,
    xTrackHeight,
    xZero,
    xMax,
    yZeroPixel,
    yMaxPixel,
    xAxisGap,
    yAxisGap,
    // GFF layout helpers passed through for redrawDotplotStage.
    yGffWidth,
    yGffSideGap,
    xGffHeight
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

function getAxisUnitForSpan(visibleSpan) {
  if (visibleSpan <= CONFIG.bpToKbThresholdBp) {
    return "bp";
  }

  if (visibleSpan <= CONFIG.kbToMbThresholdBp) {
    return "kb";
  }

  return "Mb";
}

function formatAxisValueForSpan(value, visibleSpan) {
  const unit = getAxisUnitForSpan(visibleSpan);

  if (unit === "bp") {
    return `${formatNumber(value, 0)} bp`;
  }

  if (unit === "kb") {
    return `${formatNumber(value / 1000, 1)} kb`;
  }

  return `${formatNumber(value / 1000000, 3)} Mb`;
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

  const axisY = geometry.imageHeight + 1;
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

  const axisX = geometry.imageX - 1;
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

function drawDotplotCoordinateAxes(layer, geometry, xSampleData, ySampleData) {
  drawDotplotXAxis(layer, geometry, xSampleData);
  drawDotplotYAxis(layer, geometry, ySampleData);
}

// Returns true when dotplot mode is the active viewer mode.
function isDotplotModeActive() {
  const panel = document.getElementById("dotplot-panel");
  return panel ? !panel.classList.contains("hidden") : false;
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

// Creates the dotplot Konva stage and its layers the first time dotplot mode is used.
// Also creates the persistent highlight Konva.Shape nodes and wires up all
// pointer event handlers (pointermove, pointerleave, click) for hover/pin.
// Subsequent calls are no-ops.
function initDotplotStage() {
  if (dotplotStage) {
    return;
  }
  dotplotStage = new Konva.Stage({ container: "dotplot-viewer", width: 1, height: 1 });

  dotplotImageLayer       = new Konva.Layer({ listening: false });
  dotplotTrackLayer       = new Konva.Layer({ listening: false });
  dotplotDebugLayer       = new Konva.Layer({ listening: false });
  dotplotHighlightLayer   = new Konva.Layer({ listening: false });
  dotplotInteractionLayer = new Konva.Layer();

  dotplotStage.add(dotplotImageLayer);
  dotplotStage.add(dotplotTrackLayer);
  dotplotStage.add(dotplotDebugLayer);
  dotplotStage.add(dotplotHighlightLayer);
  dotplotStage.add(dotplotInteractionLayer);

  // Persistent highlight shapes — created once, re-added to highlight layer each redraw.
  _dotplotBlockHighlightShape = new Konva.Shape({
    sceneFunc(ctx) {
      if (_dotplotBlockHighlightGeoms.length === 0) { return; }
      ctx.save();
      ctx.fillStyle = _dotplotBlockHighlightColor;
      ctx.beginPath();
      for (const r of _dotplotBlockHighlightGeoms) {
        ctx.rect(r.x, r.y, r.width, r.height);
      }
      ctx.fill();
      ctx.restore();
    },
    visible: false,
    listening: false
  });

  _dotplotSnpHighlightShape = new Konva.Shape({
    sceneFunc(ctx) {
      if (_dotplotSnpHighlightGeoms.length === 0) { return; }
      ctx.save();
      ctx.strokeStyle = _dotplotSnpHighlightColor;
      ctx.lineWidth = CONFIG.snpHighlightMinWidthPx;
      ctx.beginPath();
      for (const s of _dotplotSnpHighlightGeoms) {
        if (s.axis === "x") {
          ctx.moveTo(s.cx, s.y0);
          ctx.lineTo(s.cx, s.y1);
        } else {
          ctx.moveTo(s.x0, s.cy);
          ctx.lineTo(s.x1, s.cy);
        }
      }
      ctx.stroke();
      ctx.restore();
    },
    visible: false,
    listening: false
  });

  // Translucent blue projection bands shown on the SVG image for the
  // highlighted block's X-sample and Y-sample coordinate intervals.
  // The vertical band spans the full Y axis; the horizontal band spans the full X axis.
  _dotplotBlockIntersectionShape = new Konva.Shape({
    sceneFunc(ctx) {
      if (!_dotplotBlockIntersectionGeom) { return; }

      ctx.save();
      ctx.fillStyle = "rgba(59, 130, 246, 0.18)";

      const { vertical, horizontal } = _dotplotBlockIntersectionGeom;
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
      ctx.lineWidth = CONFIG.snpHighlightMinWidthPx;
      ctx.beginPath();

      ctx.moveTo(x, y);
      ctx.lineTo(x, yZeroPixel);

      ctx.moveTo(xZero, y);
      ctx.lineTo(x, y);

      ctx.stroke();
      ctx.restore();
    },
    visible: false,
    listening: false
  });

  // ── Pointer event handlers ─────────────────────────────────────────────────

  dotplotStage.on("pointermove", () => {
    const pointer = dotplotStage.getPointerPosition();
    if (!pointer || state.isApplyingPin) {
      return;
    }
    const resolved = resolveDotplotHoveredFeature(pointer.x, pointer.y);
    applyDotplotResolvedHover(resolved);
    // Cursor: pointer when over a feature, default otherwise.
    const container = dotplotStage.container();
    if (resolved) {
      container.style.cursor = "pointer";
    } else {
      container.style.cursor = "default";
    }
  });

  dotplotStage.on("pointerleave", () => {
    const container = dotplotStage.container();
    if (container) {
      container.style.cursor = "default";
    }
    if (state.isApplyingPin) {
      return;
    }
    applyDotplotResolvedHover(null);
  });

  dotplotStage.on("click", () => {
    if (state.isApplyingPin) {
      return;
    }
    const pointer = dotplotStage.getPointerPosition();
    if (!pointer) {
      return;
    }
    const resolved = resolveDotplotHoveredFeature(pointer.x, pointer.y);
    if (!resolved) {
      return;
    }
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

// Clears all dotplot layers and hides the stage container (used when no pair is selected).
// Also resets the hover index so it is rebuilt on next activation.
function clearDotplotStage() {
  const viewer = document.getElementById("dotplot-viewer");
  if (viewer) {
    viewer.classList.add("hidden");
  }
  _dotplotHoverIndexDirty = true;
  _lastResolvedDotplotHoverKey = null;
  if (!dotplotStage) {
    return;
  }
  dotplotImageLayer.destroyChildren();
  dotplotTrackLayer.destroyChildren();
  dotplotDebugLayer.destroyChildren();
  dotplotHighlightLayer.destroyChildren();
  dotplotStage.draw();
}

// Full batched redraw of the dotplot Konva stage.
// Uses one Konva.Shape per feature group — no one-node-per-feature.
// Visual style matches browser mode: white region, gray blocks, red SNPs, black outline.
// Also rebuilds the hover spatial index so hit-testing is always in sync with the layout.
// Computes along-axis geometry for all blocks and SNPs of one sample track,
// in a "local horizontal" coordinate system where the primary axis runs along
// the track and the cross-axis is the track height (CONFIG.trackHeight).
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
    const len    = Math.max(CONFIG.blockMinWidthPx, Math.abs(px1 - px0));
    fillRects.push({ along0, len, featureId: block.feature_id });
  }
  const snpPositions = [];
  for (const snp of snps) {
    snpPositions.push({ along: mapper(snp.pos_in_region), featureId: snp.feature_id });
  }
  return { fillRects, snpPositions };
}

function redrawDotplotStage() {
  const img = document.getElementById("dotplot-svg-img");
  if (!img || !img.complete || img.naturalWidth === 0) {
    return;
  }
  initDotplotStage();

  const geometry = computeDotplotGeometry();
  if (!geometry) {
    return;
  }

  dotplotStage.width(geometry.stageWidth);
  dotplotStage.height(geometry.stageHeight);

  const xSampleData = getSampleByName(_dotplotState.selectedX);
  const ySampleData = getSampleByName(_dotplotState.selectedY);


  // ── Image layer ────────────────────────────────────────────────────────────
  dotplotImageLayer.destroyChildren();
  dotplotImageLayer.add(new Konva.Image({
    x: geometry.imageX,
    y: geometry.imageY,
    image: img,
    width: geometry.imageWidth,
    height: geometry.imageHeight,
    listening: false
  }));

  // ── Track layer ────────────────────────────────────────────────────────────
  dotplotTrackLayer.destroyChildren();

  drawDotplotCoordinateAxes(dotplotTrackLayer, geometry, xSampleData, ySampleData);

  // Region bounds: the outer DOTPLOT_TRACK width/height = CONFIG.trackHeight + 2
  // accommodates the 1 px region border on each side (TRACK_FEATURE_INSET).
  // The inner region (xRegionH = CONFIG.trackHeight) matches the browser-mode white track rect.
  // DOTPLOT_TRACK_GAP separates the SVG image from each track region.
  const xRegionX = geometry.xZero;
  const xRegionY = geometry.imageHeight + geometry.xAxisGap + TRACK_FEATURE_INSET;
  const xRegionW = Math.max(1, geometry.xMax - geometry.xZero);
  const xRegionH = CONFIG.trackHeight;

  // Y-track region is offset right by the Y-GFF gutter + side gap so GFF tracks fit to its left.
  const yRegionX = geometry.yGffWidth + geometry.yGffSideGap + TRACK_FEATURE_INSET;
  const yRegionY = geometry.yMaxPixel;
  const yRegionW = CONFIG.trackHeight;
  const yRegionH = Math.max(1, geometry.yZeroPixel - geometry.yMaxPixel);

  // Build normalised along-axis geometry for each track using the shared helper.
  // fillRects: [{along0, len, featureId}]   — block fill positions along primary axis
  // snpPositions: [{along, featureId}]      — SNP pixel positions along primary axis
  const xGeoms = xSampleData
    ? buildTrackAlongAxisGeoms(
        xSampleData.blocks, xSampleData.snps,
        pos => mapXCoordinateToStagePx(pos, xSampleData, geometry)
      )
    : { fillRects: [], snpPositions: [] };

  const yGeoms = ySampleData
    ? buildTrackAlongAxisGeoms(
        ySampleData.blocks, ySampleData.snps,
        pos => mapYCoordinateToStagePx(pos, ySampleData, geometry)
      )
    : { fillRects: [], snpPositions: [] };

  // Stage-absolute fill/line geometry for rendering.
  // TRACK_FEATURE_INSET mirrors getFeatureY / getSnpY (1 px inset from region border).
  // X-track: horizontal — along = x, across = y.  Y-track: vertical — along = y, across = x.
  const featureH = Math.max(1, xRegionH - 2 * TRACK_FEATURE_INSET); // = CONFIG.featureHeight
  const featureW = Math.max(1, yRegionW - 2 * TRACK_FEATURE_INSET); // same for y-track

  const xBlockRects = xGeoms.fillRects.map(r => ({
    x: r.along0, y: xRegionY + TRACK_FEATURE_INSET, width: r.len, height: featureH, featureId: r.featureId
  }));
  const xSnpEntries = xGeoms.snpPositions.map(s => ({
    cx: s.along, y0: xRegionY + TRACK_FEATURE_INSET, y1: xRegionY + xRegionH - TRACK_FEATURE_INSET, featureId: s.featureId
  }));

  const yBlockRects = yGeoms.fillRects.map(r => ({
    x: yRegionX + TRACK_FEATURE_INSET, y: r.along0, width: featureW, height: r.len, featureId: r.featureId
  }));
  const ySnpEntries = yGeoms.snpPositions.map(s => ({
    cy: s.along, x0: yRegionX + TRACK_FEATURE_INSET, x1: yRegionX + yRegionW - TRACK_FEATURE_INSET, featureId: s.featureId
  }));

  // ── Build hover spatial index ───────────────────────────────────────────────
  // Store only the along-axis positions needed by resolveDotplotHoveredFeature.
  // getDotplotHighlightGeometries derives cross-axis highlight bounds from region bounds
  // + TRACK_FEATURE_INSET / TRACK_HIGHLIGHT_INSET, so they never drift from browser mode.
  const xBlocks = xBlockRects.map(r => ({ x0: r.x, x1: r.x + r.width, featureId: r.featureId }));
  xBlocks.sort((a, b) => a.x0 - b.x0);
  const xSnps = xSnpEntries.map(s => ({ cx: s.cx, featureId: s.featureId }));
  xSnps.sort((a, b) => a.cx - b.cx);

  const yBlocks = yBlockRects.map(r => ({ y0: r.y, y1: r.y + r.height, featureId: r.featureId }));
  yBlocks.sort((a, b) => a.y0 - b.y0);
  const ySnps = ySnpEntries.map(s => ({ cy: s.cy, featureId: s.featureId }));
  ySnps.sort((a, b) => a.cy - b.cy);

  _dotplotHoverIndex.xTrack = { blocks: xBlocks, snps: xSnps, regionX: xRegionX, regionY: xRegionY, regionW: xRegionW, regionH: xRegionH };
  _dotplotHoverIndex.yTrack = { blocks: yBlocks, snps: ySnps, regionX: yRegionX, regionY: yRegionY, regionW: yRegionW, regionH: yRegionH };
  _dotplotHoverIndexDirty = false;

  // ── X-track rendering ──────────────────────────────────────────────────────
  if (xSampleData) {
    // White background.
    dotplotTrackLayer.add(new Konva.Shape({
      sceneFunc(ctx) {
        ctx.fillStyle = "#ffffff";
        ctx.fillRect(xRegionX, xRegionY, xRegionW, xRegionH);
      },
      listening: false
    }));

    // Gray block fills (batched).
    if (xBlockRects.length > 0) {
      dotplotTrackLayer.add(new Konva.Shape({
        sceneFunc(ctx, shape) {
          ctx.beginPath();
          for (const r of xBlockRects) {
            ctx.rect(r.x, r.y, r.width, r.height);
          }
          ctx.fillStrokeShape(shape);
        },
        fill: CONFIG.blockFill,
        strokeWidth: 0,
        listening: false
      }));
    }

    // Red SNP lines (batched).
    if (xSnpEntries.length > 0) {
      dotplotTrackLayer.add(new Konva.Shape({
        sceneFunc(ctx, shape) {
          ctx.beginPath();
          for (const s of xSnpEntries) {
            ctx.moveTo(s.cx, s.y0);
            ctx.lineTo(s.cx, s.y1);
          }
          ctx.fillStrokeShape(shape);
        },
        stroke: CONFIG.snpColor,
        strokeWidth: CONFIG.snpMinWidthPx,
        listening: false
      }));
    }

    // Black rounded outline.
    dotplotTrackLayer.add(new Konva.Shape({
      sceneFunc(ctx, shape) {
        ctx.beginPath();
        drawRoundedRect(ctx, xRegionX, xRegionY, xRegionW, xRegionH, 2);
        ctx.fillStrokeShape(shape);
      },
      fillEnabled: false,
      stroke: "#000000",
      strokeWidth: 1,
      listening: false
    }));
  }

  // ── Y-track rendering ──────────────────────────────────────────────────────
  if (ySampleData) {
    // White background.
    dotplotTrackLayer.add(new Konva.Shape({
      sceneFunc(ctx) {
        ctx.fillStyle = "#ffffff";
        ctx.fillRect(yRegionX, yRegionY, yRegionW, yRegionH);
      },
      listening: false
    }));

    // Gray block fills (batched).
    if (yBlockRects.length > 0) {
      dotplotTrackLayer.add(new Konva.Shape({
        sceneFunc(ctx, shape) {
          ctx.beginPath();
          for (const r of yBlockRects) {
            ctx.rect(r.x, r.y, r.width, r.height);
          }
          ctx.fillStrokeShape(shape);
        },
        fill: CONFIG.blockFill,
        strokeWidth: 0,
        listening: false
      }));
    }

    // Red SNP lines (batched).
    if (ySnpEntries.length > 0) {
      dotplotTrackLayer.add(new Konva.Shape({
        sceneFunc(ctx, shape) {
          ctx.beginPath();
          for (const s of ySnpEntries) {
            ctx.moveTo(s.x0, s.cy);
            ctx.lineTo(s.x1, s.cy);
          }
          ctx.fillStrokeShape(shape);
        },
        stroke: CONFIG.snpColor,
        strokeWidth: CONFIG.snpMinWidthPx,
        listening: false
      }));
    }

    // Black rounded outline.
    dotplotTrackLayer.add(new Konva.Shape({
      sceneFunc(ctx, shape) {
        ctx.beginPath();
        drawRoundedRect(ctx, yRegionX, yRegionY, yRegionW, yRegionH, 2);
        ctx.fillStrokeShape(shape);
      },
      fillEnabled: false,
      stroke: "#000000",
      strokeWidth: 1,
      listening: false
    }));
  }

  // ── X-sample GFF tracks (horizontal, below X region) ────────────────────────
  if (xSampleData && geometry.xGffHeight > 0) {
    const xGffTracks = getSampleGffTracks(xSampleData);
    // Baseline Y: centre of each track strip, same formula as browser getGffTrackY.
    // Here panelTop equivalent = xRegionY (top of the x-track region, region height = xRegionH).
    // We place GFF tracks starting after xRegionH + GFF_TRACK.topGap below xRegionY.
    const xGffOriginY = xRegionY + xRegionH; // bottom of x-sample region (TRACK_FEATURE_INSET already counted)
    const gffRectQueuesX = new Map();

    xGffTracks.forEach((track, trackIndex) => {
      const color = getGffTrackColor(track.track_name);
      const trackY = xGffOriginY + GFF_TRACK.topGap + trackIndex * (GFF_TRACK.height + GFF_TRACK.gap);
      const baselineY = trackY + GFF_TRACK.height / 2;

      // Baseline (grey horizontal line across the full genomic range).
      dotplotTrackLayer.add(new Konva.Shape({
        sceneFunc(ctx) {
          ctx.save();
          ctx.strokeStyle = "#e5e7eb";
          ctx.lineWidth = 1;
          ctx.beginPath();
          ctx.moveTo(geometry.xZero, baselineY);
          ctx.lineTo(geometry.xMax,  baselineY);
          ctx.stroke();
          ctx.restore();
        },
        listening: false
      }));

      // Gene feature rectangles, batched by colour.
      for (const gene of track.features || []) {
        const px0 = mapXCoordinateToStagePx(gene.start_in_region, xSampleData, geometry);
        const px1 = mapXCoordinateToStagePx(gene.end_in_region,   xSampleData, geometry);
        const gx0 = Math.min(px0, px1);
        const gw  = Math.max(GFF_TRACK.minGeneWidthPx, Math.abs(px1 - px0));
        if (!gffRectQueuesX.has(color)) { gffRectQueuesX.set(color, []); }
        gffRectQueuesX.get(color).push({ x: gx0, y: trackY, width: gw, height: GFF_TRACK.height });
      }
    });

    // Flush one Konva.Shape per colour.
    for (const [color, rects] of gffRectQueuesX) {
      const rectsSnapshot = rects;
      dotplotTrackLayer.add(new Konva.Shape({
        sceneFunc(ctx, shape) {
          ctx.globalAlpha = 0.85;
          ctx.beginPath();
          for (const r of rectsSnapshot) {
            drawRoundedRect(ctx, r.x, r.y, r.width, r.height, 2);
          }
          ctx.globalAlpha = 1;
          ctx.fillStrokeShape(shape);
        },
        fill: color,
        strokeWidth: 0,
        listening: false
      }));
    }
  }

  // ── Y-sample GFF tracks (vertical, left of Y region) ────────────────────────
  if (ySampleData && geometry.yGffWidth > 0) {
    const yGffTracks = getSampleGffTracks(ySampleData);
    // X origin for track strips: they stack leftward from the y-region left edge, with a side gap.
    // yRegionX = geometry.yGffWidth + geometry.yGffSideGap + TRACK_FEATURE_INSET; strips sit to its left.
    const yGffRightEdge = geometry.yGffWidth + geometry.yGffSideGap; // left edge of y-track region (before inset)
    const gffRectQueuesY = new Map();

    yGffTracks.forEach((track, trackIndex) => {
      const color = getGffTrackColor(track.track_name);
      // Stack strips rightward from the far-left edge toward the y-region, leaving a side gap.
      // Strip 0 is nearest to the y-region.
      const stripRightX = yGffRightEdge - geometry.yGffSideGap - trackIndex * (GFF_TRACK.height + GFF_TRACK.gap);
      const trackX  = stripRightX - GFF_TRACK.height; // left edge of this strip
      const baselineX = trackX + GFF_TRACK.height / 2;

      // Baseline (grey vertical line across the genomic range).
      dotplotTrackLayer.add(new Konva.Shape({
        sceneFunc(ctx) {
          ctx.save();
          ctx.strokeStyle = "#e5e7eb";
          ctx.lineWidth = 1;
          ctx.beginPath();
          ctx.moveTo(baselineX, geometry.yMaxPixel);
          ctx.lineTo(baselineX, geometry.yZeroPixel);
          ctx.stroke();
          ctx.restore();
        },
        listening: false
      }));

      // Gene feature rectangles — rotated 90°: height=trackWidth, width=gene length.
      for (const gene of track.features || []) {
        const py0 = mapYCoordinateToStagePx(gene.start_in_region, ySampleData, geometry);
        const py1 = mapYCoordinateToStagePx(gene.end_in_region,   ySampleData, geometry);
        const gy0 = Math.min(py0, py1);
        const gh  = Math.max(GFF_TRACK.minGeneWidthPx, Math.abs(py1 - py0));
        if (!gffRectQueuesY.has(color)) { gffRectQueuesY.set(color, []); }
        gffRectQueuesY.get(color).push({ x: trackX, y: gy0, width: GFF_TRACK.height, height: gh });
      }
    });

    // Flush one Konva.Shape per colour.
    for (const [color, rects] of gffRectQueuesY) {
      const rectsSnapshot = rects;
      dotplotTrackLayer.add(new Konva.Shape({
        sceneFunc(ctx, shape) {
          ctx.globalAlpha = 0.85;
          ctx.beginPath();
          for (const r of rectsSnapshot) {
            drawRoundedRect(ctx, r.x, r.y, r.width, r.height, 2);
          }
          ctx.globalAlpha = 1;
          ctx.fillStrokeShape(shape);
        },
        fill: color,
        strokeWidth: 0,
        listening: false
      }));
    }
  }

  // ── GFF legend (bottom of stage) ───────────────────────────────────────────
  if (geometry.xGffHeight > 0) {
    const trackNames = getAllGffTrackNames();
    if (trackNames.length > 0) {
      // Legend baseline Y: bottom of the x-GFF track area.
      const legendY = geometry.stageHeight - GFF_LEGEND.height + GFF_LEGEND.topPadding;
      let legendX = geometry.xZero;
      for (const trackName of trackNames) {
        const color = getGffTrackColor(trackName);
        const textWidth = estimateTextWidth(trackName, GFF_LEGEND.fontSize);
        dotplotTrackLayer.add(new Konva.Circle({
          x: legendX + GFF_LEGEND.dotRadius,
          y: legendY + GFF_LEGEND.fontSize / 2,
          radius: GFF_LEGEND.dotRadius,
          fill: color,
          listening: false
        }));
        dotplotTrackLayer.add(new Konva.Text({
          x: legendX + GFF_LEGEND.dotRadius * 2 + GFF_LEGEND.dotTextGap,
          y: legendY,
          text: trackName,
          fontSize: GFF_LEGEND.fontSize,
          fill: "#4b5563",
          listening: false
        }));
        legendX += GFF_LEGEND.dotRadius * 2 + GFF_LEGEND.dotTextGap + textWidth + GFF_LEGEND.itemGap;
      }
    }
  }

  // ── Debug layer ────────────────────────────────────────────────────────────
  dotplotDebugLayer.destroyChildren();

  if (DOTPLOT_DEBUG_LAYOUT) {
    const { xZero, xMax, yZeroPixel, yMaxPixel, imageX, imageY, imageWidth, imageHeight } = geometry;

    // Axis boundary calibration lines — dash [10, 4], full opacity, span image area.
    dotplotDebugLayer.add(new Konva.Shape({
      sceneFunc(ctx) {
        ctx.save();
        ctx.strokeStyle = DOTPLOT_TRACK.debugColor;
        ctx.lineWidth = DOTPLOT_TRACK.debugLineWidth;
        ctx.setLineDash([10, 4]);
        ctx.beginPath();
        ctx.moveTo(xZero, imageY);                ctx.lineTo(xZero, imageY + imageHeight);
        ctx.moveTo(xMax,  imageY);                ctx.lineTo(xMax,  imageY + imageHeight);
        ctx.moveTo(imageX, yZeroPixel);           ctx.lineTo(imageX + imageWidth, yZeroPixel);
        ctx.moveTo(imageX, yMaxPixel);            ctx.lineTo(imageX + imageWidth, yMaxPixel);
        ctx.stroke();
        ctx.restore();
      },
      listening: false
    }));

    // Block boundary guide lines — dash [5, 5], half opacity, span image area.
    if (xSampleData || ySampleData) {
      const xBlockGuideXs = [];
      const yBlockGuideYs = [];

      if (xSampleData) {
        for (const block of xSampleData.blocks) {
          xBlockGuideXs.push(mapXCoordinateToStagePx(block.block_start_in_region, xSampleData, geometry));
          xBlockGuideXs.push(mapXCoordinateToStagePx(block.block_end_in_region,   xSampleData, geometry));
        }
      }
      if (ySampleData) {
        for (const block of ySampleData.blocks) {
          yBlockGuideYs.push(mapYCoordinateToStagePx(block.block_start_in_region, ySampleData, geometry));
          yBlockGuideYs.push(mapYCoordinateToStagePx(block.block_end_in_region,   ySampleData, geometry));
        }
      }

      dotplotDebugLayer.add(new Konva.Shape({
        sceneFunc(ctx) {
          ctx.save();
          ctx.strokeStyle = DOTPLOT_TRACK.debugColor;
          ctx.lineWidth = 0.8;
          ctx.globalAlpha = 0.5;
          ctx.setLineDash([5, 5]);
          ctx.beginPath();
          for (const x of xBlockGuideXs) {
            ctx.moveTo(x, imageY);
            ctx.lineTo(x, imageY + imageHeight);
          }
          for (const y of yBlockGuideYs) {
            ctx.moveTo(imageX, y);
            ctx.lineTo(imageX + imageWidth, y);
          }
          ctx.stroke();
          ctx.restore();
        },
        listening: false
      }));
    }
  }

  // ── Highlight layer ────────────────────────────────────────────────────────
  // Re-add the persistent highlight shapes and refresh their content.
  dotplotHighlightLayer.destroyChildren();
  // Intersection overlay is added first so it renders behind the track highlights.
  if (_dotplotBlockIntersectionShape) {
    dotplotHighlightLayer.add(_dotplotBlockIntersectionShape);
  }
  if (_dotplotSnpProjectionShape) {
    dotplotHighlightLayer.add(_dotplotSnpProjectionShape);
  }
  if (_dotplotBlockHighlightShape) {
    dotplotHighlightLayer.add(_dotplotBlockHighlightShape);
  }
  if (_dotplotSnpHighlightShape) {
    dotplotHighlightLayer.add(_dotplotSnpHighlightShape);
  }
  updateDotplotHighlightShapes();

  // Show the stage container now that content has been drawn.
  const viewer = document.getElementById("dotplot-viewer");
  if (viewer) {
    viewer.classList.remove("hidden");
  }
  dotplotStage.draw();

  updateDotplotAxisLabelLayout(geometry);
  // Toggle centering based on whether the stage fits the scroll container.
  // Must run after draw() so the container has its final dimensions.
  _updateDotplotScrollAlignment(geometry);
}

function updateDotplotAxisLabelLayout(geometry) {
  const xLabels = document.getElementById("dotplot-x-labels");
  const yLabels = document.getElementById("dotplot-y-labels");

  if (xLabels) {
    xLabels.style.width = `${geometry.imageWidth}px`;
    xLabels.style.marginLeft = `${geometry.imageX}px`;
  }

  if (yLabels) {
    yLabels.style.height = `${geometry.imageHeight}px`;
    yLabels.style.marginTop = `${geometry.imageY}px`;
  }
}

// Centers the dotplot stage when it fits the scroll container; left-aligns
// it when it overflows so that scrollLeft = 0 exposes the true left edge.
function _updateDotplotScrollAlignment(geometry) {
  const container = document.querySelector(".dotplot-content");
  if (!container) { return; }
  const padding = 24; // 2 × 12 px padding declared in .dotplot-content
  const available = container.clientWidth - padding;
  container.classList.toggle("dotplot-content--centered", geometry.stageWidth <= available);
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

// Resolves which feature (if any) the pointer is hovering over in dotplot mode.
// Checks x-track (horizontal strip) then y-track (vertical strip).
// SNPs are prioritized over blocks (closer distance wins; tolerance = SNP_POINTER_TOLERANCE_PX).
// Returns { featureType: "block"|"snp", featureId } or null.
function resolveDotplotHoveredFeature(pointerX, pointerY) {
  if (_dotplotHoverIndexDirty) {
    return null;
  }

  const { xTrack, yTrack } = _dotplotHoverIndex;

  // ── X-track (horizontal strip) ─────────────────────────────────────────────
  if (
    xTrack.regionH > 0 &&
    pointerX >= xTrack.regionX && pointerX <= xTrack.regionX + xTrack.regionW &&
    pointerY >= xTrack.regionY && pointerY <= xTrack.regionY + xTrack.regionH
  ) {
    // SNPs: tolerance scan around pointerX.
    const lo = lowerBoundDotplotCx(xTrack.snps, pointerX - SNP_POINTER_TOLERANCE_PX);
    let closestSnpId = null;
    let closestDist  = SNP_POINTER_TOLERANCE_PX + 1;
    for (let j = lo; j < xTrack.snps.length; j++) {
      const s = xTrack.snps[j];
      if (s.cx > pointerX + SNP_POINTER_TOLERANCE_PX) { break; }
      const d = Math.abs(pointerX - s.cx);
      if (d < closestDist) { closestDist = d; closestSnpId = s.featureId; }
    }
    if (closestSnpId !== null) {
      return { featureType: "snp", featureId: closestSnpId };
    }
    // Blocks: last block whose x0 <= pointerX, check x1.
    const bi = lowerBoundDotplotX0(xTrack.blocks, pointerX) - 1;
    if (bi >= 0 && pointerX <= xTrack.blocks[bi].x1) {
      return { featureType: "block", featureId: xTrack.blocks[bi].featureId };
    }
    return null;
  }

  // ── Y-track (vertical strip) ───────────────────────────────────────────────
  if (
    yTrack.regionH > 0 &&
    pointerX >= yTrack.regionX && pointerX <= yTrack.regionX + yTrack.regionW &&
    pointerY >= yTrack.regionY && pointerY <= yTrack.regionY + yTrack.regionH
  ) {
    // SNPs: tolerance scan around pointerY.
    const lo = lowerBoundDotplotCy(yTrack.snps, pointerY - SNP_POINTER_TOLERANCE_PX);
    let closestSnpId = null;
    let closestDist  = SNP_POINTER_TOLERANCE_PX + 1;
    for (let j = lo; j < yTrack.snps.length; j++) {
      const s = yTrack.snps[j];
      if (s.cy > pointerY + SNP_POINTER_TOLERANCE_PX) { break; }
      const d = Math.abs(pointerY - s.cy);
      if (d < closestDist) { closestDist = d; closestSnpId = s.featureId; }
    }
    if (closestSnpId !== null) {
      return { featureType: "snp", featureId: closestSnpId };
    }
    // Blocks: last block whose y0 <= pointerY, check y1.
    const bi = lowerBoundDotplotY0(yTrack.blocks, pointerY) - 1;
    if (bi >= 0 && pointerY <= yTrack.blocks[bi].y1) {
      return { featureType: "block", featureId: yTrack.blocks[bi].featureId };
    }
    return null;
  }

  return null;
}

// Applies a resolved dotplot hover, guarded by a key to avoid redundant sidebar updates.
function applyDotplotResolvedHover(resolved) {
  const key = resolved ? `${resolved.featureType}:${resolved.featureId}` : null;
  if (key === _lastResolvedDotplotHoverKey) {
    return;
  }
  _lastResolvedDotplotHoverKey = key;

  if (!resolved) {
    state.isHoveringInteractiveFeature = false;
    clearHoveredFeature();
    return;
  }

  state.isHoveringInteractiveFeature = true;
  setHoveredFeature(resolved.featureType, resolved.featureId);
}

// Computes the stage-space highlight geometry for a given feature on both dotplot tracks.
// Returns { blockGeoms: [], snpGeoms: [] } suitable for the highlight Konva.Shape nodes.
// block geoms: { x, y, width, height }
// snp geoms:   { cx, y0, y1, axis:"x" } | { cy, x0, x1, axis:"y" }
function getDotplotHighlightGeometries(featureType, featureId) {
  const blockGeoms = [];
  const snpGeoms   = [];

  if (!featureId || _dotplotHoverIndexDirty) {
    return { blockGeoms, snpGeoms };
  }

  const { xTrack, yTrack } = _dotplotHoverIndex;

  if (featureType === "block") {
    // Block highlights span the full track cross-axis with TRACK_HIGHLIGHT_INSET,
    // matching browser-mode getBlockHighlightGeometries (0.5 px inset from region border).
    for (const b of xTrack.blocks) {
      if (b.featureId === featureId) {
        blockGeoms.push({
          x:      b.x0,
          y:      xTrack.regionY + TRACK_HIGHLIGHT_INSET,
          width:  b.x1 - b.x0,
          height: xTrack.regionH - 2 * TRACK_HIGHLIGHT_INSET
        });
      }
    }
    for (const b of yTrack.blocks) {
      if (b.featureId === featureId) {
        blockGeoms.push({
          x:      yTrack.regionX + TRACK_HIGHLIGHT_INSET,
          y:      b.y0,
          width:  yTrack.regionW - 2 * TRACK_HIGHLIGHT_INSET,
          height: b.y1 - b.y0
        });
      }
    }
  } else if (featureType === "snp") {
    // SNP highlights span the same inset as feature lines (TRACK_FEATURE_INSET = 1 px).
    for (const s of xTrack.snps) {
      if (s.featureId === featureId) {
        snpGeoms.push({
          cx: s.cx,
          y0: xTrack.regionY + TRACK_FEATURE_INSET,
          y1: xTrack.regionY + xTrack.regionH - TRACK_FEATURE_INSET,
          axis: "x"
        });
      }
    }
    for (const s of yTrack.snps) {
      if (s.featureId === featureId) {
        snpGeoms.push({
          cy: s.cy,
          x0: yTrack.regionX + TRACK_FEATURE_INSET,
          x1: yTrack.regionX + yTrack.regionW - TRACK_FEATURE_INSET,
          axis: "y"
        });
      }
    }
  }

  return { blockGeoms, snpGeoms };
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
    _dotplotState.zoom = 1;

    redrawDotplotStage();

    requestAnimationFrame(() => {
      const container = document.querySelector(".dotplot-content");
      if (container) {
        container.scrollLeft = 0;
      }

      const stageEl = document.getElementById("dotplot-viewer");
      if (stageEl) {
        const stageRect = stageEl.getBoundingClientRect();
        const dotplotCenterY = window.scrollY + stageRect.top + stageRect.height / 2;
        const targetWindowScrollY = dotplotCenterY - window.innerHeight / 2;

        window.scrollTo({
          top: clampValue(
            targetWindowScrollY,
            0,
            document.documentElement.scrollHeight - window.innerHeight
          ),
          behavior: "smooth"
        });
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

  const stepPx = container.clientWidth * fraction;
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

  // Step 2 — adaptive zoom: choose the zoom that makes the rectangle occupy
  // ~50 % of the viewport on whichever axis is the binding constraint.
  // Use window dimensions for Y since we scroll the window for vertical centering.
  const TARGET_COVERAGE = 0.50;
  const viewW = container.clientWidth - 24; // subtract 2×12 px padding
  const viewH = window.innerHeight;
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

    container.scrollLeft = clampValue(
      rectCenterX - container.clientWidth / 2,
      0,
      container.scrollWidth - container.clientWidth
    );

    container.scrollTop = clampValue(
      rectCenterY - container.clientHeight / 2,
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

function _computeDotplotBlockIntersection(featureId) {
  const geometry = computeDotplotGeometry();
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

  return {
    x: x0,
    y: y0,
    width: Math.max(CONFIG.dotplotIntersectionMinSizePx, x1 - x0),
    height: Math.max(CONFIG.dotplotIntersectionMinSizePx, y1 - y0)
  };
}

// Computes the stage-space rectangle for the block projection overlay:
// the area on the SVG image that corresponds to the given block's coordinate
// interval on both the selected X and Y samples.
// Returns { x, y, width, height } in stage pixels, or null.
function _computeDotplotBlockProjection(featureId) {
  const geometry = computeDotplotGeometry();
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

  return {
    vertical: {
      x: x0,
      y: y0,
      width: Math.max(CONFIG.dotplotIntersectionMinSizePx, x1 - x0),
      height: geometry.yZeroPixel - y0
    },
    horizontal: {
      x: geometry.xZero,
      y: y0,
      width: x1 - geometry.xZero,
      height: Math.max(CONFIG.dotplotIntersectionMinSizePx, y1 - y0)
    }
  };
}


function _computeDotplotSnpProjection(featureId) {
  const geometry = computeDotplotGeometry();
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

// Updates the dotplot highlight layer to reflect the currently displayed feature
// (hover or pin), using the same color logic as browser mode.
// Must be called after updateHighlightShapes() in applyActiveDisplay(),
// and also after every full redrawDotplotStage().
function updateDotplotHighlightShapes() {
  if (!dotplotHighlightLayer || !_dotplotBlockHighlightShape || !_dotplotSnpHighlightShape) {
    return;
  }

  const displayed = getDisplayedFeature();
  const color = displayed && displayed.source === "pin"
    ? CONFIG.pinHighlightColor
    : CONFIG.hoverHighlightColor;

  let blockGeoms = [];
  let snpGeoms   = [];

  if (displayed) {
    const result = getDotplotHighlightGeometries(displayed.featureType, displayed.featureId);
    blockGeoms = result.blockGeoms;
    snpGeoms   = result.snpGeoms;
  }

  _dotplotBlockHighlightColor = color;
  _dotplotBlockHighlightGeoms = blockGeoms;
  _dotplotBlockHighlightShape.visible(blockGeoms.length > 0);

  _dotplotSnpHighlightColor = color;
  _dotplotSnpHighlightGeoms = snpGeoms;
  _dotplotSnpHighlightShape.visible(snpGeoms.length > 0);

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
    if (displayed && displayed.featureType === "snp") {
      _dotplotSnpProjectionGeom = _computeDotplotSnpProjection(displayed.featureId);
    } else {
      _dotplotSnpProjectionGeom = null;
    }
    _dotplotSnpProjectionShape.visible(_dotplotSnpProjectionGeom !== null);
  }

  dotplotHighlightLayer.batchDraw();
}

function renderDotplot() {
  const img = document.getElementById("dotplot-svg-img");

  if (!img) {
    return;
  }

  const pair = findDotplotPair(_dotplotState.selectedY, _dotplotState.selectedX);

  if (!pair) {
    updateDotplotStatusMessage("Select a compatible X/Y sample.");
    return;
  }

  _dotplotHoverIndexDirty = true;
  _lastResolvedDotplotHoverKey = null;

  img.onload = requestDotplotRedraw;
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
  _dotplotState.zoom = 1;
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
