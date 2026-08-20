const REGION_DATA = {{ REGION_DATA }};
const CONFIG = {{ CONFIG }};

Konva.pixelRatio = 1;

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

const ALIGNMENT = {
  leftMargin: 120,
  topMargin: 28,
  rowHeight: 22,
  charWidth: 14,
  minCharWidth: 6,
  maxCharWidth: 28,
  letterFontSize: 9,
  labelFontSize: 13,
  axisHeight: 24,
  bottomPadding: 28,
  scrollbarHeight: 16,
  scrollbarBottomPadding: 8,
  scrollbarMinThumbWidth: 36,
  panFraction: 0.1
};

const GFF_TRACK = {
  height: 12,
  gap: 5,
  topGap: 8,
  labelFontSize: 11,
  minGeneWidthPx: 2,
  colors: [
    "#4e79a7",
    "#f28e2b",
    "#59a14f",
    "#e15759",
    "#76b7b2",
    "#edc948",
    "#b07aa1",
    "#ff9da7",
    "#9c755f",
    "#bab0ab"
  ]
};

const GFF_LEGEND = {
  height: 22,
  dotRadius: 4,
  fontSize: 11,
  itemGap: 16,
  dotTextGap: 6,
  topPadding: 6
};

const SAMPLE_LABEL = {
  x: 24,
  rightPadding: 4,
  fontSize: 16,
  fontStyle: "bold",
  minLeftMargin: 80
};

const SNP_POINTER_TOLERANCE_PX = 8;

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

// Track dimensions: +2 (1 px inset per side) so the inner region rect matches
// CONFIG.trackHeight — the same visible height as browser-mode sample tracks.
const DOTPLOT_TRACK = {
  yTrackWidth:    CONFIG.trackHeight + 2,
  xTrackHeight:   CONFIG.trackHeight + 2,
  debugColor:     "#3b82f6",
  debugLineWidth: 1.5
};

// Track feature insets — mirror browser-mode track geometry so dotplot external
// tracks render identically to browser-mode sample tracks.
//   TRACK_FEATURE_INSET   = offset used by getFeatureY() / getSnpY() (1 px)
//   TRACK_HIGHLIGHT_INSET = offset used by getBlockHighlightGeometries() (0.5 px)
const TRACK_FEATURE_INSET   = 1;
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

const state = {
  hoveredFeatureId: null,
  hoveredFeatureType: null,
  pinnedFeatureId: null,
  pinnedFeatureType: null,
  featureGroups: new Map(),
  zoomX: 1,
  scrollX: 0,
  isDraggingViewport: false,
  dragStartPointerX: 0,
  dragStartScrollX: 0,
  isDraggingScrollbar: false,
  scrollbarDragOffsetX: 0,
  suppressHover: false,
  activeAlignmentBlockId: null,
  alignmentZoomX: 1,
  alignmentScrollX: 0,
  isDraggingAlignmentViewport: false,
  isDraggingAlignmentScrollbar: false,
  alignmentDragStartPointerX: 0,
  alignmentDragStartScrollX: 0,
  alignmentScrollbarDragOffsetX: 0,
  alignmentFocusedSnpColumn: null,
  activeKeyboardViewer: "region",
  isHoveringInteractiveFeature: false,
  isApplyingPin: false
};

const derivedData = {
  sampleOrder: [],
  allGffTrackNames: [],
  gffTrackColorByName: new Map(),
  orderedBlockFeatureIds: [],
  orderedSnpFeatureIds: [],
  kimura2pGlobalColorScaleBounds: { min: 0, max: 1 }
};

const searchIndexes = {
  blockIdToFeatureId: new Map(),
  snpKeyToFeatureId: new Map(),
  featureIdToFeatureType: new Map(),
  featureIdToRegionRange: new Map(),
  sampleByName: new Map()
};

const _searchState = {
  mode: "id",
  isOpen: false
};

const lastAlignmentRenderState = {
  blockId: undefined,
  focusedSnpColumn: undefined,
  alignmentZoomX: NaN,
  alignmentScrollX: NaN,
  containerWidth: NaN
};

const lastSidebarRenderState = {
  mode: null,
  featureId: null,
  featureType: null,
  source: null,
  isPinned: false
};

const _dotplotState = {
  selectedY: null,
  selectedX: null,
  pendingAxis: null,
  zoom: 1
};

const dotplotPairIndexes = {
  validPairsByY: new Map(),
  validPairsByX: new Map(),
  pairByYX: new Map(),
  samples: []
};

// Dotplot Konva stage and layers — initialized lazily on first dotplot activation.
let dotplotStage            = null;
let dotplotImageLayer       = null;
let dotplotTrackLayer       = null;
let dotplotDebugLayer       = null;
let dotplotHighlightLayer   = null;
let dotplotInteractionLayer = null;
let _dotplotRedrawPending   = false;

// Dotplot hover/highlight state.
// _dotplotHoverIndex holds pre-computed stage-space hit-test rectangles.
// Rebuilt only when geometry changes (pair change, image load, resize).
// Two arrays per track axis:
//   x-track blocks: [{ x0, x1, featureId }, …]  sorted by x0
//   y-track blocks: [{ y0, y1, featureId }, …]  sorted by y0
//   x-track snps:   [{ cx, featureId }, …]       sorted by cx
//   y-track snps:   [{ cy, featureId }, …]       sorted by cy
// Region bounds (regionX, regionY, regionW, regionH) are also stored and used by
// getDotplotHighlightGeometries to derive the cross-axis highlight geometry.
const _dotplotHoverIndex = {
  xTrack: { blocks: [], snps: [], regionX: 0, regionY: 0, regionW: 0, regionH: 0 },
  yTrack: { blocks: [], snps: [], regionX: 0, regionY: 0, regionW: 0, regionH: 0 }
};
let _dotplotHoverIndexDirty = true;
let _lastResolvedDotplotHoverKey = null;

// Dotplot highlight layer backing data — updated by updateDotplotHighlightShapes().
let _dotplotBlockHighlightGeoms = [];
let _dotplotBlockHighlightColor = CONFIG.hoverHighlightColor;
let _dotplotSnpHighlightGeoms   = [];
let _dotplotSnpHighlightColor   = CONFIG.hoverHighlightColor;
// These Konva.Shape nodes are created once during initDotplotStage() and
// re-added to dotplotHighlightLayer on every full redraw.
let _dotplotBlockHighlightShape    = null;
let _dotplotSnpHighlightShape      = null;
// Block intersection overlay: a single translucent blue rectangle drawn on
// the SVG image area at the intersection of the highlighted block's X- and
// Y-sample intervals.  Geometry stored as { x, y, width, height } or null.
let _dotplotBlockIntersectionShape = null;
let _dotplotBlockIntersectionGeom  = null;

let _dotplotSnpProjectionShape = null;
let _dotplotSnpProjectionGeom  = null;

function invalidateSidebarCache() {
  lastSidebarRenderState.mode = null;
  lastSidebarRenderState.featureId = null;
  lastSidebarRenderState.featureType = null;
  lastSidebarRenderState.source = null;
  lastSidebarRenderState.isPinned = false;
}

function buildFeatureGroups(data) {
  const groups = new Map();

  for (const sample of data.samples) {
    for (const block of sample.blocks) {
      const entry = {
        sample: sample.sample,
        featureType: "block",
        featureId: block.feature_id,
        info: {
          sample: sample.sample,
          block_id: block.block_id,
          block_start_in_region: block.block_start_in_region,
          block_end_in_region: block.block_end_in_region,
          block_start_in_source_seq: block.block_start_in_source_seq,
          block_end_in_source_seq: block.block_end_in_source_seq,
          length: block.block_end_in_region - block.block_start_in_region + 1
        }
      };

      if (!groups.has(block.feature_id)) {
        groups.set(block.feature_id, []);
      }
      groups.get(block.feature_id).push(entry);
    }

    for (const snp of sample.snps) {
      const entry = {
        sample: sample.sample,
        featureType: "snp",
        featureId: snp.feature_id,
        info: {
          sample: sample.sample,
          block_id: snp.block_id,
          aln_pos: snp.aln_pos,
          nt: snp.nt,
          pos_in_block: snp.pos_in_block,
          pos_in_region: snp.pos_in_region,
          pos_in_source_seq: snp.pos_in_source_seq
        }
      };

      if (!groups.has(snp.feature_id)) {
        groups.set(snp.feature_id, []);
      }
      groups.get(snp.feature_id).push(entry);
    }
  }

  return groups;
}

function getOrderedBlockFeatureIds() {
  return derivedData.orderedBlockFeatureIds;
}

function getOrderedSnpFeatureIds() {
  return derivedData.orderedSnpFeatureIds;
}

function getWrappedNeighbor(items, currentItem, direction) {
  if (items.length === 0) {
    return null;
  }

  const currentIndex = items.indexOf(currentItem);
  if (currentIndex === -1) {
    return null;
  }

  const nextIndex = (currentIndex + direction + items.length) % items.length;
  return items[nextIndex];
}

function getPinnedNavigationItems() {
  if (state.pinnedFeatureType === "block") {
    return getOrderedBlockFeatureIds();
  }

  if (state.pinnedFeatureType === "snp") {
    return getOrderedSnpFeatureIds();
  }

  return [];
}

function pinNeighborFeature(direction) {
  if (!state.pinnedFeatureId || !state.pinnedFeatureType) {
    return;
  }

  const items = getPinnedNavigationItems();
  const nextFeatureId = getWrappedNeighbor(
    items,
    state.pinnedFeatureId,
    direction
  );

  if (!nextFeatureId) {
    return;
  }

  setPinnedFeature(state.pinnedFeatureType, nextFeatureId);
}

function escapeHtml(text) {
  return String(text)
    .replaceAll("&", "&amp;")
    .replaceAll("<", "&lt;")
    .replaceAll(">", "&gt;");
}

function getSampleOrder() {
  return derivedData.sampleOrder;
}

function estimateTextWidth(text, fontSize) {
  return String(text).length * fontSize * 0.62;
}

function getLongestSampleLabelWidth() {
  return Math.max(
    ...REGION_DATA.samples.map(sample =>
      estimateTextWidth(sample.sample, SAMPLE_LABEL.fontSize)
    ),
    0
  );
}

function getLeftMargin() {
  return Math.max(
    SAMPLE_LABEL.minLeftMargin,
    SAMPLE_LABEL.x + getLongestSampleLabelWidth() + SAMPLE_LABEL.rightPadding
  );
}

function getSampleGffTracks(sample) {
  return sample.gff_tracks || [];
}

function getAllGffTrackNames() {
  return derivedData.allGffTrackNames;
}

function getGffTrackColor(trackName) {
  return derivedData.gffTrackColorByName.get(trackName) ?? "#9ca3af";
}

function getSamplePanelHeight(sample) {
  const gffTrackCount = getSampleGffTracks(sample).length;

  return CONFIG.panelHeight
    + gffTrackCount * (GFF_TRACK.height + GFF_TRACK.gap);
}

function getGffLegendHeight() {
  return getAllGffTrackNames().length > 0 ? GFF_LEGEND.height : 0;
}

function getMainViewerContentHeight() {
  const sampleHeights = REGION_DATA.samples.reduce(
    (total, sample) => total + getSamplePanelHeight(sample),
    0
  );

  return getViewerToolbarHeight()
    + CONFIG.topMargin
    + REGION_VIEWER.axisToFirstSampleGap
    + sampleHeights
    + Math.max(0, REGION_DATA.samples.length - 1) * CONFIG.panelGap
    + REGION_VIEWER.samplesToLegendGap
    + getGffLegendHeight()
    + CONFIG.bottomMargin
    + SCROLLBAR.height
    + SCROLLBAR.bottomPadding;
}

function getSamplesBottomY() {
  const sampleHeights = REGION_DATA.samples.reduce(
    (total, sample) => total + getSamplePanelHeight(sample),
    0
  );

  return getViewerToolbarHeight()
    + CONFIG.topMargin
    + REGION_VIEWER.axisToFirstSampleGap
    + sampleHeights
    + Math.max(0, REGION_DATA.samples.length - 1) * CONFIG.panelGap;
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
