
let _lastResolvedHoverKey = null;
let _hoverIndex = [];
let _hoverIndexDirty = true;
let _hoverIndexGeometryKey = "";

const SCROLLBAR = {
  height: 18,
  bottomPadding: 10,
  minThumbWidth: 36,
  trackInset: 8
};

const REGION_VIEWER = {
  axisToFirstSampleGap: 14,
  samplesToLegendGap: 10
};

const REGION_RENDERING = {
  rightMargin: 5,
  endPaddingPx: 24,
  targetTickSpacingPx: 100
};

function getHoverSpatialIndexGeometryKey() {
  return [
    getVisibleStartBp(),
    getVisibleEndBp(),
    state.zoomX,
    state.scrollX,
    getStageWidth(),
    getLeftMargin(),
    getViewerToolbarHeight()
  ].join("|");
}

function ensureHoverSpatialIndex() {
  const key = getHoverSpatialIndexGeometryKey();
  if (_hoverIndexDirty || key !== _hoverIndexGeometryKey) {
    rebuildHoverSpatialIndex();
    _hoverIndexGeometryKey = key;
    _hoverIndexDirty = false;
  }
}

function lowerBoundScreenX(arr, target) {
  let lo = 0;
  let hi = arr.length;
  while (lo < hi) {
    const mid = (lo + hi) >>> 1;
    if (arr[mid].screenX < target) {
      lo = mid + 1;
    } else {
      hi = mid;
    }
  }
  return lo;
}

function lowerBoundX0(arr, target) {
  let lo = 0;
  let hi = arr.length;
  while (lo < hi) {
    const mid = (lo + hi) >>> 1;
    if (arr[mid].x0 < target) {
      lo = mid + 1;
    } else {
      hi = mid;
    }
  }
  return lo;
}

function rebuildHoverSpatialIndex() {
  const visibleStart = getVisibleStartBp();
  const visibleEnd = getVisibleEndBp();

  _hoverIndex = REGION_DATA.samples.map((sample, i) => {
    const panelTop = computePanelTop(i);
    const trackTop = panelTop + TRACK_GEOMETRY.trackYOffset;
    const trackBottom = trackTop + TRACK_GEOMETRY.trackHeight;

    const snps = [];
    for (const snp of sample.snps) {
      if (isPositionVisible(snp.pos_in_region, visibleStart, visibleEnd)) {
        snps.push({ screenX: worldXToScreenX(snp.pos_in_region), featureId: snp.feature_id });
      }
    }
    snps.sort((a, b) => a.screenX - b.screenX);

    const blocks = [];
    for (const block of sample.blocks) {
      if (intersectsRange(block.block_start_in_region, block.block_end_in_region, visibleStart, visibleEnd)) {
        const x0 = worldXToScreenX(Math.max(block.block_start_in_region, visibleStart));
        const x1 = worldXToScreenX(Math.min(block.block_end_in_region, visibleEnd));
        const x1eff = Math.max(x0 + FEATURE_RENDERING.blockMinWidthPx, x1);
        blocks.push({ x0, x1eff, featureId: block.feature_id });
      }
    }
    blocks.sort((a, b) => a.x0 - b.x0);

    return { trackTop, trackBottom, snps, blocks };
  });
}

function resolveHoveredFeature(pointerX, pointerY) {
  ensureHoverSpatialIndex();

  for (let i = 0; i < _hoverIndex.length; i += 1) {
    const entry = _hoverIndex[i];

    if (pointerY < entry.trackTop || pointerY > entry.trackBottom) {
      continue;
    }

    const snpLo = lowerBoundScreenX(entry.snps, pointerX - SNP_POINTER_TOLERANCE_PX);
    let closestSnpFeatureId = null;
    let closestSnpDist = SNP_POINTER_TOLERANCE_PX + 1;

    for (let j = snpLo; j < entry.snps.length; j += 1) {
      const snp = entry.snps[j];
      if (snp.screenX > pointerX + SNP_POINTER_TOLERANCE_PX) {
        break;
      }
      const dist = Math.abs(pointerX - snp.screenX);
      if (dist < closestSnpDist) {
        closestSnpDist = dist;
        closestSnpFeatureId = snp.featureId;
      }
    }

    if (closestSnpFeatureId !== null) {
      return { featureType: "snp", featureId: closestSnpFeatureId };
    }

    const blockIndex = lowerBoundX0(entry.blocks, pointerX) - 1;

    if (blockIndex >= 0) {
      const block = entry.blocks[blockIndex];
      if (pointerX <= block.x1eff) {
        return { featureType: "block", featureId: block.featureId };
      }
    }

    return null;
  }

  return null;
}

function applyResolvedHover(resolved) {
  const key = resolved ? `${resolved.featureType}:${resolved.featureId}` : null;

  if (key === _lastResolvedHoverKey) {
    return;
  }

  _lastResolvedHoverKey = key;

  if (!resolved) {
    clearHoveredFeature();
    return;
  }

  setHoveredFeature(resolved.featureType, resolved.featureId);
}

function getViewerElement() {
  return document.getElementById("viewer");
}

function setBodyCursor(cursor) {
  document.body.style.cursor = cursor;
}

function setViewerCursor(cursor) {
  const viewerElement = getViewerElement();
  if (viewerElement) {
    viewerElement.style.cursor = cursor;
  }
}

function setAlignmentCursor(cursor) {
  const el = document.getElementById("alignment-viewer");
  if (el) {
    el.style.cursor = cursor;
  }
}

function getStageWidth() {
  const viewerElement = getViewerElement();
  return Math.max(1, viewerElement.clientWidth);
}

function getViewportTrackWidth() {
  return getStageWidth() - getLeftMargin() - REGION_RENDERING.rightMargin;
}

function getDrawableTrackWidth() {
  return Math.max(1, getViewportTrackWidth() - REGION_RENDERING.endPaddingPx);
}

function getVisibleBiologicalEndX() {
  return worldXToScreenX(getVisibleEndBp());
}

function getContentWidth() {
  return getDrawableTrackWidth() * state.zoomX;
}

function getMaxScrollX() {
  return Math.max(0, getContentWidth() - getDrawableTrackWidth());
}

function clampScrollX(value) {
  return Math.max(0, Math.min(getMaxScrollX(), value));
}

function normalizeScrollX() {
  state.scrollX = clampScrollX(state.scrollX);
}

function getInitialZoomX() {
  return 1;
}

function getMaxZoomX() {
  const theoretical = REGION_DATA.max_region_length / BROWSER_ZOOM.targetVisibleBp;
  return Math.min(BROWSER_ZOOM.maxZoomCap, Math.max(getInitialZoomX(), theoretical));
}

function getZoomFactor() {
  const maxZoom = getMaxZoomX();
  return Math.pow(maxZoom / getInitialZoomX(), 1 / BROWSER_ZOOM.zoomSteps);
}

function getVisibleBpSpan() {
  return REGION_DATA.max_region_length / state.zoomX;
}

function getMaxVisibleStartBp() {
  return Math.max(1, REGION_DATA.max_region_length - getVisibleBpSpan() + 1);
}

function getVisibleStartBp() {
  const maxScroll = getMaxScrollX();

  if (maxScroll <= 0) {
    return 1;
  }

  const scrollRatio = state.scrollX / maxScroll;
  return 1 + scrollRatio * (getMaxVisibleStartBp() - 1);
}

function getVisibleEndBp() {
  return Math.min(
    REGION_DATA.max_region_length,
    getVisibleStartBp() + getVisibleBpSpan() - 1
  );
}

function setVisibleStartBp(targetStartBp) {
  const clampedStart = Math.max(1, Math.min(getMaxVisibleStartBp(), targetStartBp));
  const maxScroll = getMaxScrollX();

  if (maxScroll <= 0) {
    state.scrollX = 0;
    return;
  }

  const startRatio = (clampedStart - 1) / Math.max(1, getMaxVisibleStartBp() - 1);
  state.scrollX = startRatio * maxScroll;
  normalizeScrollX();
}

function moveByViewportFraction(direction, fraction = 0.1) {
  const stepBp = getVisibleBpSpan() * fraction;
  const currentStart = getVisibleStartBp();
  setVisibleStartBp(currentStart + direction * stepBp);
  requestStageRedraw();
}

function computePanelTop(panelIndex) {
  let panelTop = getViewerToolbarHeight()
    + BROWSER_LAYOUT.topMargin
    + REGION_VIEWER.axisToFirstSampleGap;

  for (let index = 0; index < panelIndex; index += 1) {
    panelTop += getSamplePanelHeight(REGION_DATA.samples[index]) + BROWSER_LAYOUT.panelGap;
  }

  return panelTop;
}

function getMainViewerContentHeight() {
  const sampleHeights = REGION_DATA.samples.reduce(
    (total, sample) => total + getSamplePanelHeight(sample),
    0
  );

  return getViewerToolbarHeight()
    + BROWSER_LAYOUT.topMargin
    + REGION_VIEWER.axisToFirstSampleGap
    + sampleHeights
    + Math.max(0, REGION_DATA.samples.length - 1) * BROWSER_LAYOUT.panelGap
    + REGION_VIEWER.samplesToLegendGap
    + getGffLegendHeight()
    + BROWSER_LAYOUT.bottomMargin
    + SCROLLBAR.height
    + SCROLLBAR.bottomPadding;
}

function getSamplesBottomY() {
  const sampleHeights = REGION_DATA.samples.reduce(
    (total, sample) => total + getSamplePanelHeight(sample),
    0
  );

  return getViewerToolbarHeight()
    + BROWSER_LAYOUT.topMargin
    + REGION_VIEWER.axisToFirstSampleGap
    + sampleHeights
    + Math.max(0, REGION_DATA.samples.length - 1) * BROWSER_LAYOUT.panelGap;
}

function getGffLegendY() {
  return getSamplesBottomY() + REGION_VIEWER.samplesToLegendGap;
}

function drawGffTrackLegend(layer) {
  const trackNames = getAllGffTrackNames();

  if (trackNames.length === 0) {
    return;
  }

  let x = getLeftMargin();
  const y = getGffLegendY() + GFF_LEGEND.topPadding;

  for (const trackName of trackNames) {
    const color = getGffTrackColor(trackName);
    const textWidth = estimateTextWidth(trackName, GFF_LEGEND.fontSize);

    layer.add(new Konva.Circle({
      x: x + GFF_LEGEND.dotRadius,
      y: y + GFF_LEGEND.fontSize / 2,
      radius: GFF_LEGEND.dotRadius,
      fill: color,
      listening: false
    }));

    layer.add(new Konva.Text({
      x: x + GFF_LEGEND.dotRadius * 2 + GFF_LEGEND.dotTextGap,
      y,
      text: trackName,
      fontSize: GFF_LEGEND.fontSize,
      fill: "#4b5563",
      listening: false
    }));

    x += GFF_LEGEND.dotRadius * 2
      + GFF_LEGEND.dotTextGap
      + textWidth
      + GFF_LEGEND.itemGap;
  }
}

function isPointerOverSampleTrack(pointerY) {
  for (let index = 0; index < REGION_DATA.samples.length; index += 1) {
    const panelTop = computePanelTop(index);
    const trackTop = panelTop + TRACK_GEOMETRY.trackYOffset;
    const trackBottom = trackTop + TRACK_GEOMETRY.trackHeight;
    if (pointerY >= trackTop && pointerY <= trackBottom) {
      return true;
    }
  }
  return false;
}

function getScrollbarY() {
  return getMainViewerContentHeight() - SCROLLBAR.bottomPadding - SCROLLBAR.height;
}

function getWorldToScreenScale() {
  const visibleSpan = getVisibleBpSpan();
  const drawableWidth = getDrawableTrackWidth();

  if (visibleSpan <= 1) {
    return 0;
  }

  return drawableWidth / (visibleSpan - 1);
}

function worldXToScreenX(position) {
  const visibleStart = getVisibleStartBp();
  return getLeftMargin() + (position - visibleStart) * getWorldToScreenScale();
}

function formatAxisValue(value) {
  return formatGenomicCoordinate(value, getVisibleBpSpan());
}

function niceStep(value) {
  if (value <= 0) {
    return 1;
  }

  const exponent = Math.floor(Math.log10(value));
  const fraction = value / Math.pow(10, exponent);

  let niceFraction;
  if (fraction <= 1) {
    niceFraction = 1;
  } else if (fraction <= 2) {
    niceFraction = 2;
  } else if (fraction <= 5) {
    niceFraction = 5;
  } else {
    niceFraction = 10;
  }

  return niceFraction * Math.pow(10, exponent);
}

function intersectsRange(start, end, visibleStart, visibleEnd) {
  return end >= visibleStart && start <= visibleEnd;
}

function isPositionVisible(position, visibleStart, visibleEnd) {
  return position >= visibleStart && position <= visibleEnd;
}

function drawGlobalAxis(layer) {
  const x0 = getLeftMargin();
  const x1 = getVisibleBiologicalEndX();
  const axisY = getViewerToolbarHeight() + 24;

  const axis = new Konva.Line({
    points: [x0, axisY, x1, axisY],
    stroke: "#444444",
    strokeWidth: 1,
    listening: false
  });
  layer.add(axis);

  const visibleStart = getVisibleStartBp();
  const visibleEnd = getVisibleEndBp();
  const visibleSpan = Math.max(1, visibleEnd - visibleStart + 1);

  const targetPx = REGION_RENDERING.targetTickSpacingPx;
  const bpPerPixel = visibleSpan / Math.max(1, getDrawableTrackWidth());
  const rawStep = bpPerPixel * targetPx;
  const step = niceStep(rawStep);

  const firstTick = Math.ceil(visibleStart / step) * step;

  for (let value = firstTick; value <= visibleEnd; value += step) {
    const x = worldXToScreenX(value);

    const tick = new Konva.Line({
      points: [x, axisY, x, axisY + 6],
      stroke: "#444444",
      strokeWidth: 1,
      listening: false
    });

    const label = new Konva.Text({
      x: x - 30,
      y: axisY - 18,
      width: 60,
      text: formatAxisValue(value),
      fontSize: 10,
      fill: "#555555",
      align: "center",
      listening: false
    });

    layer.add(tick);
    layer.add(label);
  }
}

function drawSampleLabel(layer, panelTop, sampleName) {
  const label = new Konva.Text({
    x: SAMPLE_LABEL.x,
    y: panelTop + TRACK_GEOMETRY.trackYOffset + 4,
    width: getLeftMargin() - SAMPLE_LABEL.x - SAMPLE_LABEL.rightPadding,
    text: sampleName,
    fontSize: SAMPLE_LABEL.fontSize,
    fontStyle: SAMPLE_LABEL.fontStyle,
    fill: "#222222",
    listening: false
  });
  layer.add(label);
}

function drawSampleOutline(layer, sample, panelTop) {
  const visibleStart = getVisibleStartBp();
  const visibleEnd = getVisibleEndBp();

  if (!intersectsRange(1, sample.region_length, visibleStart, visibleEnd)) {
    return;
  }

  const clippedStart = Math.max(1, visibleStart);
  const clippedEnd = Math.min(sample.region_length, visibleEnd);

  if (clippedEnd < clippedStart) {
    return;
  }

  const x0 = worldXToScreenX(clippedStart);
  const x1 = worldXToScreenX(clippedEnd);
  const yTop = panelTop + TRACK_GEOMETRY.trackYOffset;
  const yBottom = yTop + TRACK_GEOMETRY.trackHeight;
  const isFullyVisible = clippedStart === 1 && clippedEnd === sample.region_length;

  if (isFullyVisible) {
    layer.add(new Konva.Rect({
      x: x0,
      y: yTop,
      width: Math.max(1, x1 - x0),
      height: TRACK_GEOMETRY.trackHeight,
      stroke: "black",
      strokeWidth: 1,
      fillEnabled: false,
      cornerRadius: 2,
      listening: false
    }));
    return;
  }

  layer.add(new Konva.Line({
    points: [x0, yTop, x1, yTop],
    stroke: "black",
    strokeWidth: 1,
    listening: false
  }));

  layer.add(new Konva.Line({
    points: [x0, yBottom, x1, yBottom],
    stroke: "black",
    strokeWidth: 1,
    listening: false
  }));

  if (clippedStart === 1) {
    layer.add(new Konva.Line({
      points: [x0, yTop, x0, yBottom],
      stroke: "black",
      strokeWidth: 1,
      listening: false
    }));
  }

  if (clippedEnd === sample.region_length) {
    layer.add(new Konva.Line({
      points: [x1, yTop, x1, yBottom],
      stroke: "black",
      strokeWidth: 1,
      listening: false
    }));
  }
}

function drawSampleTrackBackground(layer, sample, panelTop) {
  const visibleStart = getVisibleStartBp();
  const visibleEnd = getVisibleEndBp();

  if (!intersectsRange(1, sample.region_length, visibleStart, visibleEnd)) {
    return;
  }

  const clippedStart = Math.max(1, visibleStart);
  const clippedEnd = Math.min(sample.region_length, visibleEnd);

  if (clippedEnd < clippedStart) {
    return;
  }

  const x0 = worldXToScreenX(clippedStart);
  const x1 = worldXToScreenX(clippedEnd);
  const y = panelTop + TRACK_GEOMETRY.trackYOffset;

  layer.add(new Konva.Rect({
    x: x0,
    y,
    width: Math.max(1, x1 - x0),
    height: TRACK_GEOMETRY.trackHeight,
    fill: "#ffffff",
    strokeWidth: 0,
    listening: false
  }));
}

function getFeatureY(panelTop) {
  return panelTop + TRACK_GEOMETRY.trackYOffset + TRACK_GEOMETRY.featureInset;
}

function getSnpY(panelTop) {
  return panelTop + TRACK_GEOMETRY.trackYOffset + TRACK_GEOMETRY.featureInset;
}

function getBlockGeometry(feature, panelTop, minWidthPx) {
  const visibleStart = getVisibleStartBp();
  const visibleEnd = getVisibleEndBp();

  const clippedStart = Math.max(feature.block_start_in_region, visibleStart);
  const clippedEnd = Math.min(feature.block_end_in_region, visibleEnd);

  const x0 = worldXToScreenX(clippedStart);
  const x1 = worldXToScreenX(clippedEnd);
  const y0 = getFeatureY(panelTop);

  return {
    x: x0,
    y: y0,
    width: Math.max(minWidthPx, x1 - x0),
    height: getFeatureHeight()
  };
}

function getGffTrackY(panelTop, trackIndex) {
  return panelTop
    + TRACK_GEOMETRY.trackYOffset
    + TRACK_GEOMETRY.trackHeight
    + GFF_TRACK.topGap
    + trackIndex * (GFF_TRACK.height + GFF_TRACK.gap);
}

function drawGffTrackBaseline(layer, sample, panelTop, trackIndex) {
  const visibleStart = getVisibleStartBp();
  const visibleEnd = getVisibleEndBp();

  if (!intersectsRange(1, sample.region_length, visibleStart, visibleEnd)) {
    return;
  }

  const clippedStart = Math.max(1, visibleStart);
  const clippedEnd = Math.min(sample.region_length, visibleEnd);

  if (clippedEnd < clippedStart) {
    return;
  }

  const y = getGffTrackY(panelTop, trackIndex) + GFF_TRACK.height / 2;
  const x0 = worldXToScreenX(clippedStart);
  const x1 = worldXToScreenX(clippedEnd);

  layer.add(new Konva.Line({
    points: [x0, y, x1, y],
    stroke: "#e5e7eb",
    strokeWidth: 1,
    listening: false
  }));
}

function drawGffGeneFeature(gene, panelTop, trackIndex, color, gffRectQueues) {
  const visibleStart = getVisibleStartBp();
  const visibleEnd = getVisibleEndBp();

  if (!intersectsRange(gene.start_in_region, gene.end_in_region, visibleStart, visibleEnd)) {
    return;
  }

  const clippedStart = Math.max(gene.start_in_region, visibleStart);
  const clippedEnd = Math.min(gene.end_in_region, visibleEnd);

  if (clippedEnd < clippedStart) {
    return;
  }

  const x0 = worldXToScreenX(clippedStart);
  const x1 = worldXToScreenX(clippedEnd);
  const y = getGffTrackY(panelTop, trackIndex);

  if (!gffRectQueues.has(color)) {
    gffRectQueues.set(color, []);
  }
  gffRectQueues.get(color).push({
    x: x0,
    y,
    width: Math.max(GFF_TRACK.minGeneWidthPx, x1 - x0),
    height: GFF_TRACK.height
  });
}

function drawGffTracks(layer, sample, panelTop, gffRectQueues) {
  const tracks = getSampleGffTracks(sample);

  tracks.forEach((track, trackIndex) => {
    const color = getGffTrackColor(track.track_name);

    drawGffTrackBaseline(layer, sample, panelTop, trackIndex);

    for (const gene of track.features || []) {
      drawGffGeneFeature(gene, panelTop, trackIndex, color, gffRectQueues);
    }
  });
}

function drawSamplePanelBackground(layer, sample, panelTop) {
  const backgroundX = 4;
  const backgroundRight = getVisibleBiologicalEndX() + 8;

  layer.add(new Konva.Rect({
    x: backgroundX,
    y: panelTop - 6,
    width: Math.max(1, backgroundRight - backgroundX),
    height: getSamplePanelHeight(sample) + 4,
    fill: "#f7f8fa",
    stroke: "#e9ecef",
    strokeWidth: 1,
    cornerRadius: 8,
    listening: false
  }));
}

function drawSample(
  layer,
  outlineLayerArg,
  blockRectQueue,
  snpLineQueue,
  gffRectQueues,
  sample,
  panelIndex
) {
  const panelTop = computePanelTop(panelIndex);
  const visibleStart = getVisibleStartBp();
  const visibleEnd = getVisibleEndBp();
  drawSamplePanelBackground(layer, sample, panelTop);

  drawSampleLabel(layer, panelTop, sample.sample);
  drawSampleTrackBackground(layer, sample, panelTop);
  drawSampleOutline(outlineLayerArg, sample, panelTop);
  drawGffTracks(layer, sample, panelTop, gffRectQueues);

  for (const block of sample.blocks) {
    if (!intersectsRange(
      block.block_start_in_region,
      block.block_end_in_region,
      visibleStart,
      visibleEnd
    )) {
      continue;
    }

    blockRectQueue.push(getBlockGeometry(block, panelTop, FEATURE_RENDERING.blockMinWidthPx));
  }

  for (const snp of sample.snps) {
    if (!isPositionVisible(snp.pos_in_region, visibleStart, visibleEnd)) {
      continue;
    }

    const x = worldXToScreenX(snp.pos_in_region);
    const y0 = getSnpY(panelTop);
    snpLineQueue.push({ x, y0, y1: y0 + getSnpHeight() - 2 });
  }
}

function getScrollbarMetrics() {
  const trackX = getLeftMargin() + SCROLLBAR.trackInset;
  const trackWidth = Math.max(1, getDrawableTrackWidth() - 2 * SCROLLBAR.trackInset);
  const trackY = getScrollbarY();
  const contentWidth = getContentWidth();
  const viewportWidth = getDrawableTrackWidth();

  if (contentWidth <= viewportWidth) {
    return {
      visible: false,
      trackX,
      trackY,
      trackWidth,
      thumbX: trackX,
      thumbWidth: trackWidth
    };
  }

  const ratio = viewportWidth / contentWidth;
  const thumbWidth = Math.max(SCROLLBAR.minThumbWidth, trackWidth * ratio);
  const maxThumbTravel = Math.max(0, trackWidth - thumbWidth);
  const scrollRatio = getMaxScrollX() > 0 ? state.scrollX / getMaxScrollX() : 0;
  const thumbX = trackX + scrollRatio * maxThumbTravel;

  return {
    visible: true,
    trackX,
    trackY,
    trackWidth,
    thumbX,
    thumbWidth
  };
}

function drawScrollbar(layer) {
  const metrics = getScrollbarMetrics();

  const track = new Konva.Rect({
    x: metrics.trackX,
    y: metrics.trackY,
    width: metrics.trackWidth,
    height: SCROLLBAR.height,
    fill: "#f0f0f0",
    stroke: "#d0d0d0",
    cornerRadius: 8,
    listening: false
  });
  layer.add(track);

  const thumb = new Konva.Rect({
    x: metrics.thumbX,
    y: metrics.trackY + 1,
    width: metrics.thumbWidth,
    height: SCROLLBAR.height - 2,
    fill: metrics.visible ? "#c2c2c2" : "#e0e0e0",
    stroke: "#b4b4b4",
    cornerRadius: 7
  });

  thumb.on("mouseenter", () => {
    if (!state.isDraggingViewport && !state.isDraggingScrollbar) {
      setViewerCursor("grab");
    }
  });

  thumb.on("mouseleave", () => {
    if (!state.isDraggingViewport && !state.isDraggingScrollbar) {
      setViewerCursor("");
    }
  });

  thumb.on("pointerdown", (event) => {
    if (!metrics.visible) {
      return;
    }

    event.cancelBubble = true;
    state.isDraggingScrollbar = true;
    state.suppressHover = true;
    setBodyCursor("grabbing");

    const pointer = stage.getPointerPosition();
    state.scrollbarDragOffsetX = pointer.x - metrics.thumbX;
  });

  layer.add(thumb);

  const clickArea = new Konva.Rect({
    x: metrics.trackX,
    y: metrics.trackY,
    width: metrics.trackWidth,
    height: SCROLLBAR.height,
    fill: "rgba(0,0,0,0)"
  });

  clickArea.on("pointerdown", (event) => {
    if (!metrics.visible) {
      return;
    }

    event.cancelBubble = true;
    const pointer = stage.getPointerPosition();
    const centeredThumbX = pointer.x - metrics.thumbWidth / 2;
    setScrollFromThumbX(centeredThumbX);
    redrawStage();
  });

  layer.add(clickArea);
}

function setScrollFromThumbX(thumbX) {
  const metrics = getScrollbarMetrics();
  const maxThumbTravel = Math.max(0, metrics.trackWidth - metrics.thumbWidth);

  if (maxThumbTravel <= 0) {
    state.scrollX = 0;
    return;
  }

  const clampedThumbX = Math.max(metrics.trackX, Math.min(metrics.trackX + maxThumbTravel, thumbX));
  const thumbRatio = (clampedThumbX - metrics.trackX) / maxThumbTravel;
  state.scrollX = thumbRatio * getMaxScrollX();
  normalizeScrollX();
}

function zoomAroundViewportCenter(nextZoom) {
  const previousZoom = state.zoomX;
  const clampedZoom = Math.max(getInitialZoomX(), Math.min(getMaxZoomX(), nextZoom));

  if (clampedZoom === previousZoom) {
    return;
  }

  const oldVisibleStart = getVisibleStartBp();
  const oldVisibleSpan = getVisibleBpSpan();
  const centerBp = oldVisibleStart + oldVisibleSpan / 2;

  state.zoomX = clampedZoom;

  const newVisibleSpan = getVisibleBpSpan();
  const targetVisibleStart = centerBp - newVisibleSpan / 2;
  setVisibleStartBp(targetVisibleStart);
}

const stage = new Konva.Stage({
  container: "viewer",
  width: 1,
  height: getMainViewerContentHeight()
});

const regionLayer = new Konva.Layer({ listening: false });
const highlightLayer = new Konva.Layer({ listening: false });
const outlineLayer = new Konva.Layer({ listening: false });
const interactionLayer = new Konva.Layer();

stage.add(regionLayer);
stage.add(highlightLayer);
stage.add(outlineLayer);
stage.add(interactionLayer);

let _blockHighlightGeoms = [];
let _blockHighlightColor = FEATURE_COLORS.highlightHover;
const _blockHighlightShape = new Konva.Shape({
  sceneFunc(ctx, shape) {
    ctx.beginPath();
    for (const r of _blockHighlightGeoms) {
      ctx.rect(r.x, r.y, r.width, r.height);
    }
    ctx.fillStyle = _blockHighlightColor;
    ctx.fill();
  },
  visible: false,
  listening: false
});
highlightLayer.add(_blockHighlightShape);

let _snpHighlightGeoms = [];
let _snpHighlightColor = FEATURE_COLORS.highlightHover;
const _snpHighlightShape = new Konva.Shape({
  sceneFunc(ctx, shape) {
    ctx.beginPath();
    for (const s of _snpHighlightGeoms) {
      ctx.moveTo(s.x, s.y0);
      ctx.lineTo(s.x, s.y1);
    }
    ctx.strokeStyle = _snpHighlightColor;
    ctx.lineWidth = FEATURE_RENDERING.snpHighlightMinWidthPx;
    ctx.stroke();
  },
  visible: false,
  listening: false
});
highlightLayer.add(_snpHighlightShape);

function getBlockHighlightGeometries(featureId) {
  const visibleStart = getVisibleStartBp();
  const visibleEnd = getVisibleEndBp();
  const results = [];

  const entries = state.featureGroups.get(featureId) || [];
  for (const entry of entries) {
    const info = entry.info;
    if (!intersectsRange(info.block_start_in_region, info.block_end_in_region, visibleStart, visibleEnd)) {
      continue;
    }

    const panelIndex = viewerIndexes.sampleNameToPanelIndex.get(entry.sample);
    const panelTop = computePanelTop(panelIndex);
    const clippedStart = Math.max(info.block_start_in_region, visibleStart);
    const clippedEnd = Math.min(info.block_end_in_region, visibleEnd);
    const x0 = worldXToScreenX(clippedStart);
    const x1 = worldXToScreenX(clippedEnd);
    results.push({
      x: x0,
      y: panelTop + TRACK_GEOMETRY.trackYOffset + 0.5,
      width: Math.max(FEATURE_RENDERING.blockHighlightMinWidthPx, x1 - x0),
      height: TRACK_GEOMETRY.trackHeight - 1
    });
  }

  return results;
}

function getSnpHighlightGeometries(featureId) {
  const visibleStart = getVisibleStartBp();
  const visibleEnd = getVisibleEndBp();
  const results = [];

  const entries = state.featureGroups.get(featureId) || [];
  for (const entry of entries) {
    const info = entry.info;
    if (!isPositionVisible(info.pos_in_region, visibleStart, visibleEnd)) {
      continue;
    }

    const panelIndex = viewerIndexes.sampleNameToPanelIndex.get(entry.sample);
    const panelTop = computePanelTop(panelIndex);
    const x = worldXToScreenX(info.pos_in_region);
    const y0 = getSnpY(panelTop);
    results.push({ x, y0, y1: y0 + getSnpHeight() - 2 });
  }

  return results;
}

function updateHighlightShapes() {
  const displayed = getDisplayedFeature();
  const color = displayed && displayed.source === "pin"
    ? FEATURE_COLORS.highlightPinned
    : FEATURE_COLORS.highlightHover;

  const blockGeoms = displayed && displayed.featureType === "block"
    ? getBlockHighlightGeometries(displayed.featureId)
    : [];

  _blockHighlightColor = color;
  _blockHighlightGeoms = blockGeoms;
  _blockHighlightShape.visible(blockGeoms.length > 0);

  const snpGeoms = displayed && displayed.featureType === "snp"
    ? getSnpHighlightGeometries(displayed.featureId)
    : [];

  _snpHighlightColor = color;
  _snpHighlightGeoms = snpGeoms;
  _snpHighlightShape.visible(snpGeoms.length > 0);

  highlightLayer.batchDraw();
}

function startViewportDrag(pointerX) {
  if (getMaxScrollX() <= 0) {
    return;
  }

  state.isDraggingViewport = true;
  state.suppressHover = true;
  state.dragStartPointerX = pointerX;
  state.dragStartScrollX = state.scrollX;
  setBodyCursor("grabbing");
}

function updateViewportDrag(pointerX) {
  const deltaX = pointerX - state.dragStartPointerX;
  const worldDelta = deltaX * (getContentWidth() / getDrawableTrackWidth());
  state.scrollX = clampScrollX(state.dragStartScrollX - worldDelta);
  requestStageRedraw();
}

function updateScrollbarDrag(pointerX) {
  setScrollFromThumbX(pointerX - state.scrollbarDragOffsetX);
  requestStageRedraw();
}

function stopDrag() {
  const wasDragging = state.isDraggingViewport || state.isDraggingScrollbar;
  state.isDraggingViewport = false;
  state.isDraggingScrollbar = false;
  state.scrollbarDragOffsetX = 0;
  _lastResolvedHoverKey = null;

  if (wasDragging) {
    state.suppressHover = false;
    setBodyCursor("default");
    setViewerCursor("");
  }
}

stage.on("click", (event) => {
  if (state.isDraggingViewport || state.isDraggingScrollbar) {
    return;
  }
  const pointer = stage.getPointerPosition();
  if (!pointer) {
    return;
  }
  const scrollbarY = getScrollbarY();
  if (pointer.y >= scrollbarY) {
    return;
  }
  const resolved = resolveHoveredFeature(pointer.x, pointer.y);
  if (!resolved) {
    return;
  }
  state.isApplyingPin = true;
  state.hoveredFeatureType = null;
  state.hoveredFeatureId = null;
  _lastResolvedHoverKey = null;
  setPinnedFeature(resolved.featureType, resolved.featureId);
  requestAnimationFrame(() => {
    state.isApplyingPin = false;
  });
});

stage.on("pointerdown", (event) => {
  if (event.target !== stage) {
    return;
  }

  const pointer = stage.getPointerPosition();
  if (!pointer) {
    return;
  }

  const scrollbarY = getScrollbarY();
  if (pointer.y >= scrollbarY) {
    return;
  }

  startViewportDrag(pointer.x);
});

stage.on("pointermove", () => {
  const pointer = stage.getPointerPosition();
  if (!pointer) {
    return;
  }

  if (state.isDraggingViewport) {
    updateViewportDrag(pointer.x);
    return;
  }

  if (state.isDraggingScrollbar) {
    updateScrollbarDrag(pointer.x);
    return;
  }

  const scrollbarY = getScrollbarY();
  if (pointer.y >= scrollbarY) {
    applyResolvedHover(null);
    setViewerCursor("");
    return;
  }

  if (!state.suppressHover && !state.isApplyingPin) {
    const resolved = resolveHoveredFeature(pointer.x, pointer.y);
    applyResolvedHover(resolved);

    if (resolved) {
      setViewerCursor("pointer");
      return;
    }
  }

  if (isPointerOverSampleTrack(pointer.y)) {
    setViewerCursor("");
  } else if (getMaxScrollX() > 0) {
    setViewerCursor("grab");
  } else {
    setViewerCursor("");
  }
});

stage.on("pointerup", stopDrag);
stage.on("pointerleave", () => {
  _lastResolvedHoverKey = null;
  stopDrag();
  setViewerCursor("");
});

const alignmentStage = new Konva.Stage({
  container: "alignment-viewer",
  width: 1,
  height: 160
});

const alignmentDrawLayer = new Konva.Layer({ listening: false });
const alignmentInteractionLayer = new Konva.Layer();

alignmentStage.add(alignmentDrawLayer);
alignmentStage.add(alignmentInteractionLayer);

let _stageRedrawPending = false;
let _alignmentRedrawPending = false;
let _alignmentViewerUpdatePending = false;

function requestStageRedraw() {
  if (_stageRedrawPending) {
    return;
  }
  _stageRedrawPending = true;
  requestAnimationFrame(() => {
    _stageRedrawPending = false;
    redrawStage();
  });
}

function requestAlignmentRedraw() {
  if (_alignmentRedrawPending) {
    return;
  }
  _alignmentRedrawPending = true;
  requestAnimationFrame(() => {
    _alignmentRedrawPending = false;
    redrawAlignmentViewer();
  });
}

function requestActiveAlignmentViewerUpdate() {
  if (_alignmentViewerUpdatePending) {
    return;
  }
  _alignmentViewerUpdatePending = true;
  requestAnimationFrame(() => {
    _alignmentViewerUpdatePending = false;
    updateActiveAlignmentViewer();
  });
}

function showRenderingOverlay() {
  const overlay = document.getElementById("rendering-overlay");

  if (!overlay) {
    return;
  }

  overlay.style.display = "flex";
}

function hideRenderingOverlay() {
  const overlay = document.getElementById("rendering-overlay");
  if (overlay) {
    overlay.style.display = "none";
  }
}

function showViewerBusyOverlay(message) {
  const viewerColumn = document.getElementById("viewer-column");
  if (!viewerColumn) {
    return;
  }
  viewerColumn.style.position = "relative";
  let overlay = document.getElementById("viewer-busy-overlay");
  if (!overlay) {
    overlay = document.createElement("div");
    overlay.id = "viewer-busy-overlay";
    overlay.style.cssText = [
      "position:absolute",
      "inset:0",
      "display:flex",
      "align-items:center",
      "justify-content:center",
      "background:rgba(255,255,255,0.75)",
      "font-size:14px",
      "color:#374151",
      "z-index:100",
      "pointer-events:none",
      "user-select:none"
    ].join(";");
    viewerColumn.appendChild(overlay);
  }
  overlay.textContent = message;
  overlay.style.display = "flex";
}

function hideViewerBusyOverlay() {
  const overlay = document.getElementById("viewer-busy-overlay");
  if (overlay) {
    overlay.style.display = "none";
  }
}

function showToast(message) {
  let toast = document.getElementById("viewer-toast");

  if (!toast) {
    toast = document.createElement("div");
    toast.id = "viewer-toast";
    toast.style.cssText = [
      "position:fixed",
      "left:50%",
      "bottom:24px",
      "transform:translateX(-50%)",
      "background:rgba(17,24,39,0.92)",
      "color:white",
      "padding:10px 14px",
      "border-radius:999px",
      "font-size:13px",
      "z-index:10000",
      "opacity:0",
      "transition:opacity 180ms ease",
      "pointer-events:none"
    ].join(";");
    document.body.appendChild(toast);
  }

  toast.textContent = message;
  toast.style.opacity = "1";

  window.clearTimeout(showToast.timeoutId);
  showToast.timeoutId = window.setTimeout(() => {
    toast.style.opacity = "0";
  }, 2400);
}

function drawRoundedRect(ctx, x, y, width, height, radius) {
  const r = Math.min(radius, width / 2, height / 2);
  ctx.moveTo(x + r, y);
  ctx.lineTo(x + width - r, y);
  ctx.quadraticCurveTo(x + width, y, x + width, y + r);
  ctx.lineTo(x + width, y + height - r);
  ctx.quadraticCurveTo(x + width, y + height, x + width - r, y + height);
  ctx.lineTo(x + r, y + height);
  ctx.quadraticCurveTo(x, y + height, x, y + height - r);
  ctx.lineTo(x, y + r);
  ctx.quadraticCurveTo(x, y, x + r, y);
  ctx.closePath();
}

function redrawStage() {
  stage.width(getStageWidth());
  stage.height(getMainViewerContentHeight());

  regionLayer.destroyChildren();
  highlightLayer.destroyChildren();
  outlineLayer.destroyChildren();
  interactionLayer.destroyChildren();

  highlightLayer.add(_blockHighlightShape);
  highlightLayer.add(_snpHighlightShape);

  regionLayer.add(new Konva.Rect({
    x: 0,
    y: 0,
    width: getStageWidth(),
    height: getMainViewerContentHeight(),
    fill: "white",
    listening: false
  }));

  drawGlobalAxis(regionLayer);

  const blockRectQueue = [];
  const snpLineQueue = [];
  const gffRectQueues = new Map();

  REGION_DATA.samples.forEach((sample, index) => {
    drawSample(
      regionLayer,
      outlineLayer,
      blockRectQueue,
      snpLineQueue,
      gffRectQueues,
      sample,
      index
    );
  });

  if (blockRectQueue.length > 0) {
    regionLayer.add(new Konva.Shape({
      sceneFunc(ctx, shape) {
        ctx.beginPath();
        for (const rect of blockRectQueue) {
          ctx.rect(rect.x, rect.y, rect.width, rect.height);
        }
        ctx.fillStrokeShape(shape);
      },
      fill: FEATURE_COLORS.block,
      strokeWidth: 0,
      listening: false
    }));
  }

  if (snpLineQueue.length > 0) {
    regionLayer.add(new Konva.Shape({
      sceneFunc(ctx, shape) {
        ctx.beginPath();
        for (const seg of snpLineQueue) {
          ctx.moveTo(seg.x, seg.y0);
          ctx.lineTo(seg.x, seg.y1);
        }
        ctx.fillStrokeShape(shape);
      },
      stroke: FEATURE_COLORS.snp,
      strokeWidth: FEATURE_RENDERING.snpMinWidthPx,
      listening: false
    }));
  }

  gffRectQueues.forEach((rects, color) => {
    regionLayer.add(new Konva.Shape({
      sceneFunc(ctx, shape) {
        ctx.beginPath();
        for (const r of rects) {
          drawRoundedRect(ctx, r.x, r.y, r.width, r.height, 2);
        }
        ctx.fillStrokeShape(shape);
      },
      fill: color,
      opacity: 0.85,
      strokeWidth: 0,
      listening: false
    }));
  });

  drawGffTrackLegend(regionLayer);
  drawScrollbar(interactionLayer);
  reapplyDisplayIfVisible();
  stage.draw();
  ensureHoverSpatialIndex();
}
