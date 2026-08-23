const REGION_DATA = {{ REGION_DATA }};

// Konva.pixelRatio = 1;

const FEATURE_COLORS = {
  block: "#d9d9d9",
  snp: "#d62728",
  rejectedSnp: "#666666",
  highlightPinned: "rgb(0,120,255)",
  highlightHover: "#93c5fd"
};

const TRACK_GEOMETRY = {
  trackHeight: 28,
  trackYOffset: 12,
  featureInset: 1
};

const FEATURE_RENDERING = {
  blockMinWidthPx: 1,
  blockHighlightMinWidthPx: 2,
  snpMinWidthPx: 1,
  snpHighlightMinWidthPx: 2
};

const COORDINATE_FORMAT = {
  bpToKbThresholdBp: 10_000,
  kbToMbThresholdBp: 1_000_000
};

const BROWSER_LAYOUT = {
  topMargin: 30,
  bottomMargin: 10,
  panelGap: 10,
  panelBottomSpace: 12
};

const BROWSER_ZOOM = {
  targetVisibleBp: 2_000,
  maxZoomCap: 10_000,
  zoomSteps: 16
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

const GFF_HOVER = {
  maxVisibleBp: 1_000_000,
  maxVisibleGenesPerSample: 100
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

function formatNumber(value, decimals = 0) {
  const fixed = value.toFixed(decimals);
  const parts = fixed.split(".");

  if (parts.length === 1) {
    return parts[0];
  }

  const trimmedFraction = parts[1].replace(/0+$/, "");
  if (trimmedFraction === "") {
    return parts[0];
  }

  return `${parts[0]}.${trimmedFraction}`;
}

function getGenomicCoordinateUnit(referenceValue) {
  if (referenceValue <= COORDINATE_FORMAT.bpToKbThresholdBp) {
    return "bp";
  }

  if (referenceValue <= COORDINATE_FORMAT.kbToMbThresholdBp) {
    return "kb";
  }

  return "Mb";
}

function formatGenomicCoordinate(value, referenceValue) {
  const unit = getGenomicCoordinateUnit(referenceValue);

  if (unit === "bp") {
    return `${formatNumber(value, 0)} bp`;
  }

  if (unit === "kb") {
    return `${formatNumber(value / 1000, 1)} kb`;
  }

  return `${formatNumber(value / 1000000, 3)} Mb`;
}

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
  isApplyingPin: false,
  showRejectedSnps: true
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

const viewerIndexes = {
  sampleNameToPanelIndex: new Map(),
  blockSnpByBlockId: new Map()
};

let _hasWarnedMissingSnpResult = false;

function getSnpResult(featureId) {
  const result = REGION_DATA.snp_results?.[featureId];

  if (!result && !_hasWarnedMissingSnpResult) {
    console.warn(
      `Missing workflow result for SNP feature "${featureId}"; keeping it visible.`
    );
    _hasWarnedMissingSnpResult = true;
  }

  return result || null;
}

function isRejectedSnp(featureId) {
  return getSnpResult(featureId)?.final_status === "FAIL";
}

function shouldDisplaySnp(featureId) {
  return state.showRejectedSnps || !isRejectedSnp(featureId);
}

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

// Dotplot rendering surfaces — initialized lazily on first dotplot activation.
// The matrix and both genomic context panes deliberately use separate stages so
// browser-native scrolling can keep the context panes frozen without Konva work.
let dotplotStage             = null;
let dotplotImageLayer        = null;
let dotplotDebugLayer        = null;
let dotplotHighlightLayer    = null;
let dotplotXStage            = null;
let dotplotXTrackLayer       = null;
let dotplotXHighlightLayer   = null;
let dotplotYStage            = null;
let dotplotYTrackLayer       = null;
let dotplotYHighlightLayer   = null;
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

// Dotplot context highlight backing data — updated by updateDotplotHighlightShapes().
let _dotplotXBlockHighlightGeoms = [];
let _dotplotXBlockHighlightColor = FEATURE_COLORS.highlightHover;
let _dotplotXSnpHighlightGeoms   = [];
let _dotplotXSnpHighlightColor   = FEATURE_COLORS.highlightHover;
let _dotplotYBlockHighlightGeoms = [];
let _dotplotYBlockHighlightColor = FEATURE_COLORS.highlightHover;
let _dotplotYSnpHighlightGeoms   = [];
let _dotplotYSnpHighlightColor   = FEATURE_COLORS.highlightHover;
// These Konva.Shape nodes are created once and re-added to their corresponding
// context highlight layer on every full redraw.
let _dotplotXBlockHighlightShape = null;
let _dotplotXSnpHighlightShape   = null;
let _dotplotYBlockHighlightShape = null;
let _dotplotYSnpHighlightShape   = null;
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

function buildViewerIndexes() {
  const featureGroups = new Map();
  const blockPositions = new Map();
  const snpPositions = new Map();
  const trackNameSet = new Set();
  const blockSnpAccumulators = new Map();

  searchIndexes.blockIdToFeatureId.clear();
  searchIndexes.snpKeyToFeatureId.clear();
  searchIndexes.featureIdToFeatureType.clear();
  searchIndexes.featureIdToRegionRange.clear();
  searchIndexes.sampleByName.clear();
  viewerIndexes.sampleNameToPanelIndex.clear();
  viewerIndexes.blockSnpByBlockId.clear();

  derivedData.sampleOrder = [];

  for (const [panelIndex, sample] of REGION_DATA.samples.entries()) {
    const sampleName = sample.sample;
    derivedData.sampleOrder.push(sampleName);
    searchIndexes.sampleByName.set(sampleName, sample);
    viewerIndexes.sampleNameToPanelIndex.set(sampleName, panelIndex);

    for (const track of (sample.gff_tracks || [])) {
      trackNameSet.add(track.track_name);
    }

    for (const block of sample.blocks) {
      const entry = {
        sample: sampleName,
        featureType: "block",
        featureId: block.feature_id,
        info: {
          sample: sampleName,
          block_id: block.block_id,
          block_start_in_region: block.block_start_in_region,
          block_end_in_region: block.block_end_in_region,
          block_start_in_source_seq: block.block_start_in_source_seq,
          block_end_in_source_seq: block.block_end_in_source_seq,
          length: block.block_end_in_region - block.block_start_in_region + 1
        }
      };

      if (!featureGroups.has(block.feature_id)) {
        featureGroups.set(block.feature_id, []);
      }
      featureGroups.get(block.feature_id).push(entry);

      const numericBlockId = Number(block.block_id);
      if (!searchIndexes.blockIdToFeatureId.has(numericBlockId)) {
        searchIndexes.blockIdToFeatureId.set(numericBlockId, block.feature_id);
      }
      searchIndexes.featureIdToFeatureType.set(block.feature_id, "block");
      if (!searchIndexes.featureIdToRegionRange.has(block.feature_id)) {
        searchIndexes.featureIdToRegionRange.set(block.feature_id, {
          start: block.block_start_in_region,
          end: block.block_end_in_region
        });
      }
      if (!blockPositions.has(block.feature_id)) {
        blockPositions.set(block.feature_id, Number(block.block_start_in_region));
      }
    }

    for (const snp of sample.snps) {
      const entry = {
        sample: sampleName,
        featureType: "snp",
        featureId: snp.feature_id,
        info: {
          sample: sampleName,
          block_id: snp.block_id,
          aln_pos: snp.aln_pos,
          nt: snp.nt,
          pos_in_block: snp.pos_in_block,
          pos_in_region: snp.pos_in_region,
          pos_in_source_seq: snp.pos_in_source_seq
        }
      };

      if (!featureGroups.has(snp.feature_id)) {
        featureGroups.set(snp.feature_id, []);
      }
      featureGroups.get(snp.feature_id).push(entry);

      const snpKey = `${snp.block_id}:${snp.aln_pos}`;
      if (!searchIndexes.snpKeyToFeatureId.has(snpKey)) {
        searchIndexes.snpKeyToFeatureId.set(snpKey, snp.feature_id);
      }
      searchIndexes.featureIdToFeatureType.set(snp.feature_id, "snp");
      if (!searchIndexes.featureIdToRegionRange.has(snp.feature_id)) {
        searchIndexes.featureIdToRegionRange.set(snp.feature_id, {
          start: snp.pos_in_region,
          end: snp.pos_in_region
        });
      }
      if (!snpPositions.has(snp.feature_id)) {
        snpPositions.set(snp.feature_id, Number(snp.pos_in_region));
      }

      const alignmentPosition = Number(snp.aln_pos);
      if (Number.isNaN(alignmentPosition) || alignmentPosition < 1) {
        continue;
      }

      const blockId = String(snp.block_id);
      if (!blockSnpAccumulators.has(blockId)) {
        blockSnpAccumulators.set(blockId, {
          alignmentColumns: new Set(),
          navigationItemsByFeatureId: new Map()
        });
      }

      const blockSnpAccumulator = blockSnpAccumulators.get(blockId);
      const columnIndex = alignmentPosition - 1;
      blockSnpAccumulator.alignmentColumns.add(columnIndex);
      if (!blockSnpAccumulator.navigationItemsByFeatureId.has(snp.feature_id)) {
        blockSnpAccumulator.navigationItemsByFeatureId.set(snp.feature_id, {
          featureId: snp.feature_id,
          columnIndex
        });
      }
    }
  }

  state.featureGroups = featureGroups;
  derivedData.allGffTrackNames = [...trackNameSet].sort();
  derivedData.gffTrackColorByName = new Map(
    derivedData.allGffTrackNames.map(function(name, index) {
      return [name, GFF_TRACK.colors[index % GFF_TRACK.colors.length]];
    })
  );
  derivedData.orderedBlockFeatureIds = [...blockPositions.entries()]
    .sort(function(a, b) { return a[1] - b[1]; })
    .map(function(entry) { return entry[0]; });
  derivedData.orderedSnpFeatureIds = [...snpPositions.entries()]
    .sort(function(a, b) { return a[1] - b[1]; })
    .map(function(entry) { return entry[0]; });

  for (const [blockId, accumulator] of blockSnpAccumulators) {
    viewerIndexes.blockSnpByBlockId.set(blockId, {
      alignmentColumns: accumulator.alignmentColumns,
      navigationItems: [...accumulator.navigationItemsByFeatureId.values()]
        .sort((left, right) => left.columnIndex - right.columnIndex)
    });
  }

  initializeKimura2pGlobalColorScaleBounds();
}

function initializeKimura2pGlobalColorScaleBounds() {
  const k2pValues = [];
  Object.entries(REGION_DATA.kimura2p_matrices || {}).forEach(function([blockId, matrix]) {
    if (!matrix || !matrix.values) {
      return;
    }
    matrix.values.forEach(function(row, rowIndex) {
      row.forEach(function(value, colIndex) {
        if (colIndex <= rowIndex) {
          return;
        }
        const numericValue = Number(value);
        if (!Number.isNaN(numericValue)) {
          k2pValues.push(numericValue);
        }
      });
    });
  });

  if (k2pValues.length === 0) {
    console.warn("Kimura 2P color scale: no numeric off-diagonal values found.");
    derivedData.kimura2pGlobalColorScaleBounds = { min: 0, max: 1 };
  } else {
    derivedData.kimura2pGlobalColorScaleBounds = {
      min: Math.min(...k2pValues),
      max: Math.max(...k2pValues)
    };
  }
}

function getOrderedBlockFeatureIds() {
  return derivedData.orderedBlockFeatureIds;
}

function getOrderedSnpFeatureIds() {
  return derivedData.orderedSnpFeatureIds.filter(shouldDisplaySnp);
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

function getFeatureHeight() {
  return TRACK_GEOMETRY.trackHeight - 2 * TRACK_GEOMETRY.featureInset;
}

function getSnpHeight() {
  return TRACK_GEOMETRY.trackHeight;
}

function getBaseSamplePanelHeight() {
  return TRACK_GEOMETRY.trackYOffset
    + TRACK_GEOMETRY.trackHeight
    + BROWSER_LAYOUT.panelBottomSpace;
}

function getSamplePanelHeight(sample) {
  const gffTrackCount = getSampleGffTracks(sample).length;

  return getBaseSamplePanelHeight()
    + gffTrackCount * (GFF_TRACK.height + GFF_TRACK.gap);
}

function getGffLegendHeight() {
  return getAllGffTrackNames().length > 0 ? GFF_LEGEND.height : 0;
}
