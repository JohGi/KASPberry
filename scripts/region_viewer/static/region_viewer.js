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

function formatDistanceValue(value) {
  if (value === null || value === undefined || Number.isNaN(Number(value))) {
    return "NA";
  }

  const numericValue = Number(value);
  const normalizedValue = Math.abs(numericValue) < 1e-9 ? 0 : numericValue;

  if (Math.abs(normalizedValue) >= 10) {
    return formatNumber(normalizedValue, 1);
  }

  if (Math.abs(normalizedValue) >= 1) {
    return formatNumber(normalizedValue, 2);
  }

  return formatNumber(normalizedValue, 4);
}

function getMatrixColorScaleBounds(matrix) {
  const values = [];

  matrix.values.forEach((row, rowIndex) => {
    row.forEach((value, colIndex) => {
      if (rowIndex === colIndex) {
        return;
      }

      const numericValue = Number(value);
      if (!Number.isNaN(numericValue)) {
        values.push(numericValue);
      }
    });
  });

  if (values.length === 0) {
    return { min: 0, max: 1 };
  }

  return {
    min: Math.min(...values),
    max: Math.max(...values)
  };
}

function getKimura2pGlobalColorScaleBounds() {
  return derivedData.kimura2pGlobalColorScaleBounds;
}

function interpolateChannel(start, end, ratio) {
  return Math.round(start + (end - start) * ratio);
}

function interpolateRgb(start, end, ratio) {
  return [
    interpolateChannel(start[0], end[0], ratio),
    interpolateChannel(start[1], end[1], ratio),
    interpolateChannel(start[2], end[2], ratio)
  ];
}

function matrixCellColor(value, minValue, maxValue, isDiagonal) {
  if (isDiagonal) {
    return "rgb(102, 190, 125)";
  }

  const numericValue = Number(value);
  if (Number.isNaN(numericValue) || maxValue <= minValue) {
    return "#ffffff";
  }

  const ratio = Math.max(0, Math.min(1, (numericValue - minValue) / (maxValue - minValue)));

  const green = [102, 190, 125];
  const yellow = [255, 235, 132];
  const orange = [248, 172, 89];

  const rgb = ratio <= 0.5
    ? interpolateRgb(green, yellow, ratio / 0.5)
    : interpolateRgb(yellow, orange, (ratio - 0.5) / 0.5);

  return `rgb(${rgb[0]}, ${rgb[1]}, ${rgb[2]})`;
}

function getDistanceMatrixTooltip(matrix) {
  if (matrix.source === "kimura2p") {
    return `Values are substitutions per base (Kimura 2-parameter).\nComputed from unmasked block alignments using EMBOSS distmat (-nucmethod 2), scaled from per 100 bases.`;
  }

  return "";
}

function renderInfoTooltip(text) {
  if (!text) {
    return "";
  }

  return `<span class="info-tooltip" data-tooltip="${escapeHtml(text)}">Info</span>`;
}

const VIEWER_MODE_INFO_TOOLTIPS = {
  browser: [
    "This view shows locally collinear blocks retained across samples in gray, and SNPs retained after filtering in red.",
    "",
    "Locally collinear blocks were first identified across the input samples with SibeliaZ, then filtered to retain blocks that were present exactly once in each sample, had consistent strand orientation, and passed the configured minimum block length.",
    "",
    "For each retained block, sequences were extracted, hard-masked with RepeatMasker, and aligned with MAFFT. Biallelic SNPs were then called from the masked block alignments with SeqTUI, using the configured minimum flank length, i.e., the minimum number of perfectly aligned bases required on both sides of each SNP, with no indel, additional SNP, or N. Depending on the analysis configuration, the displayed SNPs may also have been further restricted to a marker subset or to group-discriminant variants.",
    "",
    "References: SibeliaZ, Minkin & Medvedev 2020; RepeatMasker, Tarailo-Graovac & Chen 2009; MAFFT, Katoh & Standley 2013; SeqTUI, Ranwez 2026."
  ].join("\n"),

  dotplot: [
    "This view shows pairwise dotplots between available sample pairs.",
    "",
    "For each selected pair, the corresponding sequences were first aligned with Minimap2 (-x asm5), then visualized as a pairwise dotplot with blastn2dotplots.",
    "",
    "References: Minimap2, Li 2018; blastn2dotplots, Okuno et al. 2025."
  ].join("\n")
};


function updateViewerModeInfoTooltip(mode) {
  const infoTooltip = document.getElementById("viewer-mode-info-tooltip");

  if (!infoTooltip) {
    return;
  }

  infoTooltip.dataset.tooltip = VIEWER_MODE_INFO_TOOLTIPS[mode] || "";
}

function renderDistanceMatrix(matrix) {
  if (!matrix || !matrix.labels || !matrix.values) {
    return `
      <div class="distance-matrix-card">
        <p class="hint">No distance matrix available.</p>
      </div>
    `;
  }

  const labels = matrix.labels;
  const tooltip = getDistanceMatrixTooltip(matrix);
  const colorScale = matrix.source === "kimura2p"
    ? getKimura2pGlobalColorScaleBounds()
    : getMatrixColorScaleBounds(matrix);

  let html = `
    <div class="distance-matrix-card">
      <p class="distance-matrix-title">
        ${escapeHtml(matrix.title || "Distance matrix")}
        ${renderInfoTooltip(tooltip)}
      </p>
      <div class="distance-matrix-scroll">
        <table class="distance-matrix-table">
          <thead>
            <tr>
              <th></th>
  `;

  for (const label of labels) {
    html += `<th>${escapeHtml(label)}</th>`;
  }

  html += `
            </tr>
          </thead>
          <tbody>
  `;

  labels.forEach((rowLabel, rowIndex) => {
    html += `<tr><th class="row-label">${escapeHtml(rowLabel)}</th>`;

    labels.forEach((_colLabel, colIndex) => {
      if (colIndex < rowIndex) {
        html += `<td class="matrix-empty"></td>`;
        return;
      }

      const value = matrix.values[rowIndex]?.[colIndex];
      const backgroundColor = matrixCellColor(
        value,
        colorScale.min,
        colorScale.max,
        rowIndex === colIndex
      );

      html += `<td style="background-color: ${backgroundColor};">${escapeHtml(formatDistanceValue(value))}</td>`;
    });

    html += `</tr>`;
  });

  html += `
          </tbody>
        </table>
      </div>
    </div>
  `;

  return html;
}

function renderGlobalSummaryStats() {
  const globalStats = REGION_DATA.summary_stats?.global;

  if (!globalStats) {
    return `
      <div class="summary-card">
        <h3>Global statistics</h3>
        <p class="hint">No summary statistics available.</p>
      </div>
    `;
  }

  const entries = [
    ["Kept blocks", globalStats.n_blocks_kept],
    ["Smallest block (bp)", globalStats.min_block_len_bp],
    ["Largest block (bp)", globalStats.max_block_len_bp],
    ["Mean block length (bp)", globalStats.mean_block_len_bp],
    ["Kept SNPs", globalStats.n_snps_kept]
  ];

  let html = `
    <div class="summary-card">
      <h3>Global statistics</h3>
      <div class="kv">
  `;

  for (const [label, value] of entries) {
    html += `<div class="key">${escapeHtml(label)}</div><div>${escapeHtml(value)}</div>`;
  }

  html += `
      </div>
    </div>
  `;

  return html;
}

function formatSnpGroups(settings) {
  const groups = settings.snp_groups || {};

  const entries = Object.entries(groups);

  if (entries.length === 0) {
    return "none";
  }

  return entries
    .map(([group, samples]) => `${group}: ${samples.join(", ")}`)
    .join(" | ");
}
// function formatSnpGroups(settings) {
//   const groupA = settings.snp_group_a || [];
//   const groupB = settings.snp_group_b || [];

//   if (groupA.length === 0 && groupB.length === 0) {
//     return "none";
//   }

//   if (groupA.length === 0) {
//     return `all other samples vs ${groupB.join(", ")}`;
//   }

//   if (groupB.length === 0) {
//     return `${groupA.join(", ")} vs all other samples`;
//   }

//   return `${groupA.join(", ")} vs ${groupB.join(", ")}`;
// }

function renderAnalysisSettings() {
  const settings = REGION_DATA.analysis_settings;
  const container = document.getElementById("analysis-settings-content");

  if (!container) {
    return;
  }

  if (!settings || Object.keys(settings).length === 0) {
    return;
  }

  const minBlock = settings.minimum_block_length_bp !== null && settings.minimum_block_length_bp !== undefined
    ? `${settings.minimum_block_length_bp} bp`
    : "NA";

  const minFlank = settings.minimum_snp_flank_bp !== null && settings.minimum_snp_flank_bp !== undefined
    ? `${settings.minimum_snp_flank_bp} bp`
    : "NA";

  const snpGroups = formatSnpGroups(settings);

  const entries = [
    ["Minimum block length", minBlock],
    ["Minimum SNP flank", minFlank],
    ["SNP filtering groups", snpGroups]
  ];

  let html = `<div class="summary-card"><div class="kv">`;

  for (const [label, value] of entries) {
    html += `<div class="key">${escapeHtml(label)}</div><div>${escapeHtml(value)}</div>`;
  }

  html += `</div></div>`;

  container.innerHTML = html;
}

function renderSampleRegionStats() {
  const sampleStats = REGION_DATA.summary_stats?.samples;

  if (!sampleStats || Object.keys(sampleStats).length === 0) {
    return `
      <div class="summary-card">
        <h3>Sample region statistics</h3>
        <p class="hint">No per-sample region statistics available.</p>
      </div>
    `;
  }

  let html = `<div class="summary-card"><h3>Sample region statistics</h3>`;

  for (const sampleName of getSampleOrder()) {
    const stats = sampleStats[sampleName];
    html += `<div class="sample-card"><h3>${escapeHtml(sampleName)}</h3>`;

    if (!stats) {
      html += `<p class="hint">No data.</p>`;
    } else {
      const entries = [
        ["Region length (bp)", stats.region_length_bp],
        ["Covered by blocks (%)", `${formatNumber(Number(stats.covered_pct_of_region), 2)}%`]
      ];

      html += `<div class="kv">`;
      for (const [label, value] of entries) {
        html += `<div class="key">${escapeHtml(label)}</div><div>${escapeHtml(String(value))}</div>`;
      }
      html += `</div>`;
    }

    html += `</div>`;
  }

  html += `</div>`;
  return html;
}

function getMaskedNStats(blockId, sampleName) {
  return REGION_DATA.masked_block_n_stats?.[String(blockId)]?.[sampleName] || null;
}

function formatFeatureInfoEntries(featureType, info) {
  if (featureType === "block") {
    return [
      ["Coords in region", `${info.block_start_in_region}-${info.block_end_in_region}`],
      [
        "Coords in source seq",
        `${info.block_start_in_source_seq}-${info.block_end_in_source_seq}`
      ],
      ["Length", String(info.length)]
    ];
  }

  return [
    ["Allele", String(info.nt)],
    ["Pos in region", String(info.pos_in_region)],
    ["Pos in source seq", String(info.pos_in_source_seq)]
  ];
}

function getBlockHighlightGeometries(featureId) {
  const visibleStart = getVisibleStartBp();
  const visibleEnd = getVisibleEndBp();
  const results = [];
  for (let i = 0; i < REGION_DATA.samples.length; i += 1) {
    const sample = REGION_DATA.samples[i];
    for (const block of sample.blocks) {
      if (block.feature_id !== featureId) {
        continue;
      }
      if (!intersectsRange(block.block_start_in_region, block.block_end_in_region, visibleStart, visibleEnd)) {
        continue;
      }
      const panelTop = computePanelTop(i);
      const clippedStart = Math.max(block.block_start_in_region, visibleStart);
      const clippedEnd = Math.min(block.block_end_in_region, visibleEnd);
      const x0 = worldXToScreenX(clippedStart);
      const x1 = worldXToScreenX(clippedEnd);
      results.push({
        x: x0,
        y: panelTop + CONFIG.trackY + 0.5,
        width: Math.max(CONFIG.blockHighlightMinWidthPx, x1 - x0),
        height: CONFIG.trackHeight - 1
      });
    }
  }
  return results;
}

function getSnpHighlightGeometries(featureId) {
  const visibleStart = getVisibleStartBp();
  const visibleEnd = getVisibleEndBp();
  const results = [];
  for (let i = 0; i < REGION_DATA.samples.length; i += 1) {
    const sample = REGION_DATA.samples[i];
    for (const snp of sample.snps) {
      if (snp.feature_id !== featureId) {
        continue;
      }
      if (!isPositionVisible(snp.pos_in_region, visibleStart, visibleEnd)) {
        continue;
      }
      const panelTop = computePanelTop(i);
      const x = worldXToScreenX(snp.pos_in_region);
      const y0 = getSnpY(panelTop);
      results.push({ x, y0, y1: y0 + CONFIG.snpHeight - 2 });
    }
  }
  return results;
}

function updateHighlightShapes() {
  const displayed = getDisplayedFeature();
  const color = displayed && displayed.source === "pin"
    ? CONFIG.pinHighlightColor
    : CONFIG.hoverHighlightColor;

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

function getDisplayedFeature() {
  if (state.hoveredFeatureId && state.hoveredFeatureType) {
    return {
      featureId: state.hoveredFeatureId,
      featureType: state.hoveredFeatureType,
      source: "hover"
    };
  }

  if (state.pinnedFeatureId && state.pinnedFeatureType) {
    return {
      featureId: state.pinnedFeatureId,
      featureType: state.pinnedFeatureType,
      source: "pin"
    };
  }

  return null;
}

function updateAnalysisSettingsVisibility() {
  const analysisPanel = document.getElementById("analysis-settings-sidebar");

  if (!analysisPanel) {
    return;
  }

  const hasDisplayedFeature = Boolean(getDisplayedFeature());
  analysisPanel.classList.toggle("hidden", hasDisplayedFeature);
  syncSidebarHeightToViewerColumn();
}

function applyActiveDisplay() {
  updateFeatureNavigationButtons();
  updateHighlightShapes();
  if (isDotplotModeActive()) {
    updateDotplotHighlightShapes();
  }
  updateAnalysisSettingsVisibility();

  const displayed = getDisplayedFeature();

  if (!displayed) {
    renderSidebarDefault();
    requestActiveAlignmentViewerUpdate();
    stage.batchDraw();
    return;
  }

  renderFeatureSidebar(
    displayed.featureType,
    displayed.featureId,
    displayed.source === "pin"
  );
  requestActiveAlignmentViewerUpdate();
  stage.batchDraw();
}

function renderSidebarDefault() {
  if (lastSidebarRenderState.mode === "default") {
    return;
  }

  lastSidebarRenderState.mode = "default";
  lastSidebarRenderState.featureId = null;
  lastSidebarRenderState.featureType = null;
  lastSidebarRenderState.source = null;
  lastSidebarRenderState.isPinned = false;

  const sidebar = document.getElementById("sidebar");
  sidebar.innerHTML = `
    <div class="sidebar-header">
      <h2>Region overview</h2>
    </div>
    <div class="sidebar-section">
      ${renderDistanceMatrix(REGION_DATA.mash_matrix)}
    </div>
    <div class="sidebar-section">
      ${renderGlobalSummaryStats()}
    </div>
    <div class="sidebar-section">
      ${renderSampleRegionStats()}
    </div>
  `;
}

function renderSidebarHeader(title, isPinned) {
  return `
    <div class="sidebar-header">
      <div>
        <h2>${escapeHtml(title)}</h2>
        ${isPinned ? '<div class="pin-badge">Pinned</div>' : ""}
      </div>
      ${state.pinnedFeatureId ? '<button id="sidebar-unpin" class="sidebar-close" type="button" aria-label="Unpin feature">✕</button>' : ""}
    </div>
  `;
}

function attachSidebarUnpinHandler() {
  const unpinButton = document.getElementById("sidebar-unpin");
  if (unpinButton) {
    unpinButton.addEventListener("click", () => {
      clearPinnedFeature();
    });
  }
}

function updateFeatureNavigationButtons() {
  const hasPinnedFeature = Boolean(state.pinnedFeatureId && state.pinnedFeatureType);

  const navGroup = document.getElementById("feature-nav-group");
  const prevBtn = document.getElementById("feature-prev");
  const nextBtn = document.getElementById("feature-next");
  const centerBtn = document.getElementById("feature-center");
  const zoomGroup = document.getElementById("viewer-zoom-group");

  if (navGroup) {
    navGroup.classList.toggle("hidden", !hasPinnedFeature);
  }
  if (prevBtn) {
    prevBtn.classList.toggle("hidden", !hasPinnedFeature);
  }
  if (nextBtn) {
    nextBtn.classList.toggle("hidden", !hasPinnedFeature);
  }
  if (centerBtn) {
    centerBtn.classList.toggle("hidden", !hasPinnedFeature);
  }
  if (zoomGroup) {
    zoomGroup.classList.toggle("has-separator", hasPinnedFeature);
  }
}

function renderBlockSidebar(featureId, isPinned) {
  const sidebar = document.getElementById("sidebar");
  const entries = state.featureGroups.get(featureId) || [];

  if (entries.length === 0) {
    renderSidebarDefault();
    return;
  }

  const firstInfo = entries[0].info;
  const blockId = String(firstInfo.block_id);
  const matrix = REGION_DATA.kimura2p_matrices?.[blockId];

  let html = `
    ${renderSidebarHeader("Block", isPinned)}
    <p class="hint"><b>ID:</b> ${escapeHtml(blockId)}</p>
    <div class="sidebar-section">
      ${renderDistanceMatrix(matrix)}
    </div>
  `;

  for (const sampleName of getSampleOrder()) {
    const entry = entries.find(item => item.sample === sampleName);
    html += `<div class="sample-card"><h3>${escapeHtml(sampleName)}</h3>`;

    if (!entry) {
      html += `<p class="hint">No corresponding feature in this sample.</p>`;
    } else {
      const nStats = getMaskedNStats(blockId, sampleName);
      const formattedEntries = formatFeatureInfoEntries("block", entry.info);

      if (nStats) {
        formattedEntries.push(["Total N (%)", `${formatNumber(Number(nStats.masked_n_pct), 2)}%`]);
        formattedEntries.push(["Repeat/TE N (%)", `${formatNumber(Number(nStats.repeat_masked_n_pct), 2)}%`]);
      } else {
        formattedEntries.push(["Total N (%)", "NA"]);
        formattedEntries.push(["Repeat/TE N (%)", "NA"]);
      }

      html += '<div class="kv">';

      for (const [label, value] of formattedEntries) {
        html += `<div class="key">${escapeHtml(label)}</div><div>${escapeHtml(value)}</div>`;
      }

      html += "</div>";
    }

    html += "</div>";
  }

  sidebar.innerHTML = html;
  attachSidebarUnpinHandler();
}

function renderSnpSidebar(featureId, isPinned) {
  const sidebar = document.getElementById("sidebar");
  const entries = state.featureGroups.get(featureId) || [];

  if (entries.length === 0) {
    renderSidebarDefault();
    return;
  }

  const firstInfo = entries[0].info;
  const title = `${firstInfo.block_id}:${firstInfo.aln_pos}`;

  let html = `${renderSidebarHeader("SNP", isPinned)}<p class="hint"><b>ID:</b> ${escapeHtml(title)}</p>`;

  for (const sampleName of getSampleOrder()) {
    const entry = entries.find(item => item.sample === sampleName);
    html += `<div class="sample-card"><h3>${escapeHtml(sampleName)}</h3>`;

    if (!entry) {
      html += `<p class="hint">No corresponding feature in this sample.</p>`;
    } else {
      html += '<div class="kv">';

      const formattedEntries = formatFeatureInfoEntries("snp", entry.info);
      for (const [label, value] of formattedEntries) {
        html += `<div class="key">${escapeHtml(label)}</div><div>${escapeHtml(value)}</div>`;
      }

      html += "</div>";
    }

    html += "</div>";
  }

  sidebar.innerHTML = html;
  attachSidebarUnpinHandler();
}

function renderFeatureSidebar(featureType, featureId, isPinned) {
  const source = isPinned ? "pin" : "hover";

  if (
    lastSidebarRenderState.mode === "feature" &&
    lastSidebarRenderState.featureId === featureId &&
    lastSidebarRenderState.featureType === featureType &&
    lastSidebarRenderState.source === source &&
    lastSidebarRenderState.isPinned === isPinned
  ) {
    return;
  }

  lastSidebarRenderState.mode = "feature";
  lastSidebarRenderState.featureId = featureId;
  lastSidebarRenderState.featureType = featureType;
  lastSidebarRenderState.source = source;
  lastSidebarRenderState.isPinned = isPinned;

  if (featureType === "block") {
    renderBlockSidebar(featureId, isPinned);
    return;
  }

  renderSnpSidebar(featureId, isPinned);
}

function setHoveredFeature(featureType, featureId) {
  state.hoveredFeatureType = featureType;
  state.hoveredFeatureId = featureId;
  applyActiveDisplay();
}

function clearHoveredFeature() {
  state.hoveredFeatureType = null;
  state.hoveredFeatureId = null;
  applyActiveDisplay();
}

function setPinnedFeature(featureType, featureId) {
  state.pinnedFeatureType = featureType;
  state.pinnedFeatureId = featureId;
  applyActiveDisplay();
}

function clearPinnedFeature() {
  state.pinnedFeatureType = null;
  state.pinnedFeatureId = null;
  invalidateSidebarCache();
  applyActiveDisplay();
}

function reapplyDisplayIfVisible() {
  applyActiveDisplay();
}

let _lastResolvedHoverKey = null;
let _hoverIndex = [];
let _hoverIndexDirty = true;
let _hoverIndexGeometryKey = "";

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
    const trackTop = panelTop + CONFIG.trackY;
    const trackBottom = trackTop + CONFIG.trackHeight;

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
        const x1eff = Math.max(x0 + CONFIG.blockMinWidthPx, x1);
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
    state.isHoveringInteractiveFeature = false;
    clearHoveredFeature();
    return;
  }

  state.isHoveringInteractiveFeature = true;
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

function resetViewerCursors() {
  setViewerCursor("");
  setAlignmentCursor("");
}

function getStageWidth() {
  const viewerElement = getViewerElement();
  return Math.max(1, viewerElement.clientWidth);
}

function getViewportTrackWidth() {
  return getStageWidth() - getLeftMargin() - CONFIG.rightMargin;
}

function getDrawableTrackWidth() {
  return Math.max(1, getViewportTrackWidth() - CONFIG.endPaddingPx);
}

function getVisibleBiologicalEndX() {
  return worldXToScreenX(getVisibleEndBp());
}

function getVisibleBiologicalWidth() {
  return Math.max(1, getVisibleBiologicalEndX() - getLeftMargin());
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
  const theoretical = REGION_DATA.max_region_length / CONFIG.targetVisibleBp;
  return Math.min(CONFIG.maxZoomCap, Math.max(getInitialZoomX(), theoretical));
}

function getZoomFactor() {
  const maxZoom = getMaxZoomX();
  return Math.pow(maxZoom / getInitialZoomX(), 1 / CONFIG.zoomSteps);
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
    + CONFIG.topMargin
    + REGION_VIEWER.axisToFirstSampleGap;

  for (let index = 0; index < panelIndex; index += 1) {
    panelTop += getSamplePanelHeight(REGION_DATA.samples[index]) + CONFIG.panelGap;
  }

  return panelTop;
}

function isPointerOverSampleTrack(pointerY) {
  for (let index = 0; index < REGION_DATA.samples.length; index += 1) {
    const panelTop = computePanelTop(index);
    const trackTop = panelTop + CONFIG.trackY;
    const trackBottom = trackTop + CONFIG.trackHeight;
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

function getAxisUnit() {
  const visibleSpan = getVisibleBpSpan();

  if (visibleSpan <= CONFIG.bpToKbThresholdBp) {
    return "bp";
  }

  if (visibleSpan <= CONFIG.kbToMbThresholdBp) {
    return "kb";
  }

  return "Mb";
}

function formatAxisValue(value) {
  const unit = getAxisUnit();

  if (unit === "bp") {
    return `${formatNumber(value, 0)} bp`;
  }

  if (unit === "kb") {
    return `${formatNumber(value / 1000, 1)} kb`;
  }

  return `${formatNumber(value / 1000000, 3)} Mb`;
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

  const targetPx = CONFIG.targetTickSpacingPx;
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
    y: panelTop + CONFIG.trackY + 4,
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
  const yTop = panelTop + CONFIG.trackY;
  const yBottom = yTop + CONFIG.trackHeight;
  const isFullyVisible = clippedStart === 1 && clippedEnd === sample.region_length;

  if (isFullyVisible) {
    layer.add(new Konva.Rect({
      x: x0,
      y: yTop,
      width: Math.max(1, x1 - x0),
      height: CONFIG.trackHeight,
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
  const y = panelTop + CONFIG.trackY;

  layer.add(new Konva.Rect({
    x: x0,
    y,
    width: Math.max(1, x1 - x0),
    height: CONFIG.trackHeight,
    fill: "#ffffff",
    strokeWidth: 0,
    listening: false
  }));
}

function getFeatureY(panelTop) {
  return panelTop + CONFIG.trackY + 1;
}

function getSnpY(panelTop) {
  return panelTop + CONFIG.trackY + 1;
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
    height: CONFIG.featureHeight
  };
}

function getGffTrackY(panelTop, trackIndex) {
  return panelTop
    + CONFIG.trackY
    + CONFIG.trackHeight
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

    blockRectQueue.push(getBlockGeometry(block, panelTop, CONFIG.blockMinWidthPx));
  }

  for (const snp of sample.snps) {
    if (!isPositionVisible(snp.pos_in_region, visibleStart, visibleEnd)) {
      continue;
    }

    const x = worldXToScreenX(snp.pos_in_region);
    const y0 = getSnpY(panelTop);
    snpLineQueue.push({ x, y0, y1: y0 + CONFIG.snpHeight - 2 });
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
let _blockHighlightColor = CONFIG.hoverHighlightColor;
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
let _snpHighlightColor = CONFIG.hoverHighlightColor;
const _snpHighlightShape = new Konva.Shape({
  sceneFunc(ctx, shape) {
    ctx.beginPath();
    for (const s of _snpHighlightGeoms) {
      ctx.moveTo(s.x, s.y0);
      ctx.lineTo(s.x, s.y1);
    }
    ctx.strokeStyle = _snpHighlightColor;
    ctx.lineWidth = CONFIG.snpHighlightMinWidthPx;
    ctx.stroke();
  },
  visible: false,
  listening: false
});
highlightLayer.add(_snpHighlightShape);

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
      fill: CONFIG.blockFill,
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
      stroke: CONFIG.snpColor,
      strokeWidth: CONFIG.snpMinWidthPx,
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

function getAlignmentContainer() {
  return document.getElementById("alignment-viewer");
}

function getAlignmentSubtitle() {
  return document.getElementById("alignment-subtitle");
}

function getAlignmentPanel() {
  return document.getElementById("alignment-panel");
}

function getBlockIdFromFeature(featureType, featureId) {
  if (featureType !== "block") {
    return null;
  }

  const entries = state.featureGroups.get(featureId) || [];
  if (entries.length === 0) {
    return null;
  }

  return String(entries[0].info.block_id);
}

function getFirstFeatureInfo(featureId) {
  const entries = state.featureGroups.get(featureId) || [];

  if (entries.length === 0) {
    return null;
  }

  return entries[0].info;
}

function getActiveAlignmentTarget() {
  if (state.hoveredFeatureId && state.hoveredFeatureType === "block") {
    return {
      blockId: getBlockIdFromFeature("block", state.hoveredFeatureId),
      focusColumn: null
    };
  }

  if (state.pinnedFeatureId && state.pinnedFeatureType === "snp") {
    const info = getFirstFeatureInfo(state.pinnedFeatureId);

    if (!info) {
      return null;
    }

    return {
      blockId: String(info.block_id),
      focusColumn: Number(info.aln_pos) - 1
    };
  }

  if (state.pinnedFeatureId && state.pinnedFeatureType === "block") {
    return {
      blockId: getBlockIdFromFeature("block", state.pinnedFeatureId),
      focusColumn: null
    };
  }

  return null;
}

function getAlignmentData(blockId) {
  if (!blockId) {
    return null;
  }

  return REGION_DATA.block_alignments?.[String(blockId)] || null;
}

function getAlignmentSampleNames(alignmentData) {
  if (!alignmentData) {
    return [];
  }

  const sampleOrder = getSampleOrder();
  const availableSamples = new Set(Object.keys(alignmentData));
  return sampleOrder.filter(sampleName => availableSamples.has(sampleName));
}

function getAlignmentLength(alignmentData) {
  if (!alignmentData) {
    return 0;
  }

  const sequences = Object.values(alignmentData);
  if (sequences.length === 0) {
    return 0;
  }

  return sequences[0].length;
}

function getAlignmentStageWidth() {
  const container = getAlignmentContainer();
  return Math.max(1, container.clientWidth);
}

function getAlignmentStageHeight(sampleCount) {
  return ALIGNMENT.topMargin
    + ALIGNMENT.axisHeight
    + sampleCount * ALIGNMENT.rowHeight
    + ALIGNMENT.bottomPadding
    + ALIGNMENT.scrollbarHeight
    + ALIGNMENT.scrollbarBottomPadding;
}

function getAlignmentViewportWidth() {
  return Math.max(1, getAlignmentStageWidth() - ALIGNMENT.leftMargin - 12);
}

function getAlignmentCharWidth() {
  return Math.max(
    ALIGNMENT.minCharWidth,
    Math.min(ALIGNMENT.maxCharWidth, ALIGNMENT.charWidth * state.alignmentZoomX)
  );
}

function getAlignmentContentWidth(alignmentLength) {
  return alignmentLength * getAlignmentCharWidth();
}

function getAlignmentMaxScrollX(alignmentLength) {
  return Math.max(
    0,
    getAlignmentContentWidth(alignmentLength) - getAlignmentViewportWidth()
  );
}

function clampAlignmentScrollX(value, alignmentLength) {
  return Math.max(0, Math.min(getAlignmentMaxScrollX(alignmentLength), value));
}

function normalizeAlignmentScrollX(alignmentLength) {
  state.alignmentScrollX = clampAlignmentScrollX(
    state.alignmentScrollX,
    alignmentLength
  );
}

function centerAlignmentOnColumn(columnIndex, alignmentLength) {
  if (
    columnIndex === null
    || columnIndex === undefined
    || Number.isNaN(Number(columnIndex))
    || alignmentLength === 0
  ) {
    return;
  }

  const charWidth = getAlignmentCharWidth();
  const targetCenterX = columnIndex * charWidth + charWidth / 2;

  state.alignmentScrollX = clampAlignmentScrollX(
    targetCenterX - getAlignmentViewportWidth() / 2,
    alignmentLength
  );
}

function getVisibleAlignmentColumnRange(alignmentLength) {
  const charWidth = getAlignmentCharWidth();
  const firstColumn = Math.max(0, Math.floor(state.alignmentScrollX / charWidth));
  const visibleColumnCount = Math.ceil(getAlignmentViewportWidth() / charWidth) + 2;
  const lastColumn = Math.min(alignmentLength, firstColumn + visibleColumnCount);

  return {
    firstColumn,
    lastColumn
  };
}

function alignmentColumnToX(columnIndex) {
  return ALIGNMENT.leftMargin
    + columnIndex * getAlignmentCharWidth()
    - state.alignmentScrollX;
}

function alignmentColumnCenterToX(columnIndex) {
  return alignmentColumnToX(columnIndex) + getAlignmentCharWidth() / 2;
}

function getBaseFill(base) {
  const normalizedBase = String(base).toUpperCase();

  if (normalizedBase === "A") {
    return "#7fc97f";
  }

  if (normalizedBase === "C") {
    return "#80b1d3";
  }

  if (normalizedBase === "G") {
    return "#fdb462";
  }

  if (normalizedBase === "T") {
    return "#fb8072";
  }

  if (normalizedBase === "-") {
    return "#e5e7eb";
  }

  if (normalizedBase === "N") {
    return "#d1d5db";
  }

  return "#f3f4f6";
}

function getBaseTextFill(base) {
  if (String(base) === "-") {
    return "#6b7280";
  }

  return "#111827";
}

function formatAlignmentAxisValue(value) {
  if (value <= CONFIG.bpToKbThresholdBp) {
    return `${formatNumber(value, 0)} bp`;
  }

  if (value <= CONFIG.kbToMbThresholdBp) {
    return `${formatNumber(value / 1000, 1)} kb`;
  }

  return `${formatNumber(value / 1000000, 3)} Mb`;
}

function drawAlignmentAxis(layer, alignmentLength, snpColumns) {
  const y = ALIGNMENT.topMargin;
  const x0 = ALIGNMENT.leftMargin;
  const x1 = ALIGNMENT.leftMargin + getAlignmentViewportWidth();

  layer.add(new Konva.Line({
    points: [x0, y, x1, y],
    stroke: "#444444",
    strokeWidth: 1,
    listening: false
  }));

  const charWidth = getAlignmentCharWidth();
  const bpPerPixel = 1 / Math.max(1, charWidth);
  const rawStep = bpPerPixel * 90;
  const step = niceStep(rawStep);
  const range = getVisibleAlignmentColumnRange(alignmentLength);
  const firstTick = Math.ceil((range.firstColumn + 1) / step) * step;

  for (let value = firstTick; value <= range.lastColumn; value += step) {
    const x = alignmentColumnCenterToX(value - 1);

    layer.add(new Konva.Line({
      points: [x, y, x, y + 6],
      stroke: "#444444",
      strokeWidth: 1,
      listening: false
    }));

    layer.add(new Konva.Text({
      x: x - 34,
      y: y - 18,
      width: 68,
      text: formatAlignmentAxisValue(value),
      fontSize: 10,
      fill: "#555555",
      align: "center",
      listening: false
    }));
  }
  const visibleSnpColumns = [...snpColumns].filter(
    columnIndex => columnIndex >= range.firstColumn && columnIndex < range.lastColumn
  );

  for (const columnIndex of visibleSnpColumns) {
    const x = alignmentColumnCenterToX(columnIndex);

    layer.add(new Konva.Text({
      x: x - 6,
      y: y + 12,
      width: 12,
      text: "*",
      fontSize: 13,
      fontStyle: "bold",
      fill: "#111827",
      align: "center",
      listening: false
    }));
  }
}

function getBlockSnpAlignmentColumns(blockId) {
  const snpColumns = new Set();

  if (!blockId) {
    return snpColumns;
  }

  for (const sample of REGION_DATA.samples) {
    for (const snp of sample.snps) {
      if (String(snp.block_id) !== String(blockId)) {
        continue;
      }

      const alnPos = Number(snp.aln_pos);
      if (!Number.isNaN(alnPos) && alnPos >= 1) {
        snpColumns.add(alnPos - 1);
      }
    }
  }

  return snpColumns;
}

function getBlockSnpColumns(blockId) {
  const columns = [];

  if (!blockId) {
    return columns;
  }

  for (const sample of REGION_DATA.samples) {
    for (const snp of sample.snps) {
      if (String(snp.block_id) !== String(blockId)) {
        continue;
      }

      const alnPos = Number(snp.aln_pos);
      if (!Number.isNaN(alnPos) && alnPos >= 1) {
        columns.push(alnPos - 1);
      }
    }
  }

  return [...new Set(columns)].sort((left, right) => left - right);
}

function getBlockSnpNavigationItems(blockId) {
  const itemsByFeatureId = new Map();

  if (!blockId) {
    return [];
  }

  for (const sample of REGION_DATA.samples) {
    for (const snp of sample.snps) {
      if (String(snp.block_id) !== String(blockId)) {
        continue;
      }

      const alnPos = Number(snp.aln_pos);
      if (Number.isNaN(alnPos) || alnPos < 1) {
        continue;
      }

      if (!itemsByFeatureId.has(snp.feature_id)) {
        itemsByFeatureId.set(snp.feature_id, {
          featureId: snp.feature_id,
          columnIndex: alnPos - 1
        });
      }
    }
  }

  return [...itemsByFeatureId.values()]
    .sort((left, right) => left.columnIndex - right.columnIndex);
}

function getCurrentAlignmentCenterColumn() {
  return (state.alignmentScrollX + getAlignmentViewportWidth() / 2)
    / getAlignmentCharWidth();
}

function getNearestSnpColumnIndex(snpColumns) {
  if (snpColumns.length === 0) {
    return -1;
  }

  const currentColumn = state.alignmentFocusedSnpColumn !== null
    ? state.alignmentFocusedSnpColumn
    : getCurrentAlignmentCenterColumn();

  let bestIndex = 0;
  let bestDistance = Math.abs(snpColumns[0] - currentColumn);

  for (let index = 1; index < snpColumns.length; index += 1) {
    const distance = Math.abs(snpColumns[index] - currentColumn);

    if (distance < bestDistance) {
      bestDistance = distance;
      bestIndex = index;
    }
  }

  return bestIndex;
}

function focusNeighborAlignmentSnp(direction) {
  const blockId = state.activeAlignmentBlockId;
  const alignmentData = getAlignmentData(blockId);
  const alignmentLength = getAlignmentLength(alignmentData);
  const snpItems = getBlockSnpNavigationItems(blockId);

  if (alignmentLength === 0 || snpItems.length === 0) {
    return;
  }

  const snpColumns = snpItems.map(item => item.columnIndex);
  const currentIndex = getNearestSnpColumnIndex(snpColumns);

  if (currentIndex === -1) {
    return;
  }

  const nextIndex = (currentIndex + direction + snpItems.length) % snpItems.length;
  const nextItem = snpItems[nextIndex];

  state.alignmentFocusedSnpColumn = nextItem.columnIndex;
  setPinnedFeature("snp", nextItem.featureId);
}

function drawAlignmentRows(layer, alignmentData, sampleNames, alignmentLength, snpColumns) {
  const range = getVisibleAlignmentColumnRange(alignmentLength);
  const charWidth = getAlignmentCharWidth();
  const baseY = ALIGNMENT.topMargin + ALIGNMENT.axisHeight;

  sampleNames.forEach((sampleName, rowIndex) => {
    const rowY = baseY + rowIndex * ALIGNMENT.rowHeight;
    const sequence = alignmentData[sampleName] || "";

    layer.add(new Konva.Text({
      x: 10,
      y: rowY + 4,
      width: ALIGNMENT.leftMargin - 18,
      text: sampleName,
      fontSize: ALIGNMENT.labelFontSize,
      fontStyle: "bold",
      fill: "#222222",
      listening: false
    }));

    for (let columnIndex = range.firstColumn; columnIndex < range.lastColumn; columnIndex += 1) {
      const base = sequence[columnIndex] || "";
      const x = alignmentColumnToX(columnIndex);
      const isSnpColumn = snpColumns.has(columnIndex);

      layer.add(new Konva.Rect({
        x,
        y: rowY,
        width: Math.max(1, charWidth),
        height: ALIGNMENT.rowHeight - 1,
        fill: getBaseFill(base),
        stroke: "#ffffff",
        strokeWidth: 0.5,
        listening: false
      }));

      if (charWidth >= 9) {
        layer.add(new Konva.Text({
          x,
          y: rowY,
          width: charWidth,
          height: ALIGNMENT.rowHeight - 1,
          text: base,
          fontSize: isSnpColumn
            ? ALIGNMENT.letterFontSize + 2
            : ALIGNMENT.letterFontSize,
          fontStyle: isSnpColumn ? "bold" : "normal",
          fill: getBaseTextFill(base),
          align: "center",
          verticalAlign: "middle",
          listening: false
        }));
      }
    }
  });
}

function getAlignmentScrollbarY(sampleCount) {
  return getAlignmentStageHeight(sampleCount)
    - ALIGNMENT.scrollbarBottomPadding
    - ALIGNMENT.scrollbarHeight;
}

function getAlignmentScrollbarMetrics(alignmentLength, sampleCount) {
  const trackX = ALIGNMENT.leftMargin;
  const trackY = getAlignmentScrollbarY(sampleCount);
  const trackWidth = getAlignmentViewportWidth();
  const contentWidth = getAlignmentContentWidth(alignmentLength);
  const viewportWidth = getAlignmentViewportWidth();

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
  const thumbWidth = Math.max(ALIGNMENT.scrollbarMinThumbWidth, trackWidth * ratio);
  const maxThumbTravel = Math.max(0, trackWidth - thumbWidth);
  const maxScroll = getAlignmentMaxScrollX(alignmentLength);
  const scrollRatio = maxScroll > 0 ? state.alignmentScrollX / maxScroll : 0;
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

function setAlignmentScrollFromThumbX(thumbX, alignmentLength, sampleCount) {
  const metrics = getAlignmentScrollbarMetrics(alignmentLength, sampleCount);
  const maxThumbTravel = Math.max(0, metrics.trackWidth - metrics.thumbWidth);

  if (maxThumbTravel <= 0) {
    state.alignmentScrollX = 0;
    return;
  }

  const clampedThumbX = Math.max(
    metrics.trackX,
    Math.min(metrics.trackX + maxThumbTravel, thumbX)
  );
  const thumbRatio = (clampedThumbX - metrics.trackX) / maxThumbTravel;
  state.alignmentScrollX = thumbRatio * getAlignmentMaxScrollX(alignmentLength);
  normalizeAlignmentScrollX(alignmentLength);
}

function drawAlignmentScrollbar(layer, alignmentLength, sampleCount) {
  const metrics = getAlignmentScrollbarMetrics(alignmentLength, sampleCount);

  layer.add(new Konva.Rect({
    x: metrics.trackX,
    y: metrics.trackY,
    width: metrics.trackWidth,
    height: ALIGNMENT.scrollbarHeight,
    fill: "#f0f0f0",
    stroke: "#d0d0d0",
    cornerRadius: 8,
    listening: false
  }));

  const thumb = new Konva.Rect({
    x: metrics.thumbX,
    y: metrics.trackY + 1,
    width: metrics.thumbWidth,
    height: ALIGNMENT.scrollbarHeight - 2,
    fill: metrics.visible ? "#c2c2c2" : "#e0e0e0",
    stroke: "#b4b4b4",
    cornerRadius: 7
  });

  thumb.on("pointerdown", (event) => {
    if (!metrics.visible) {
      return;
    }

    event.cancelBubble = true;
    state.isDraggingAlignmentScrollbar = true;
    setBodyCursor("grabbing");

    const pointer = alignmentStage.getPointerPosition();
    state.alignmentScrollbarDragOffsetX = pointer.x - metrics.thumbX;
  });

  layer.add(thumb);

  const clickArea = new Konva.Rect({
    x: metrics.trackX,
    y: metrics.trackY,
    width: metrics.trackWidth,
    height: ALIGNMENT.scrollbarHeight,
    fill: "rgba(0,0,0,0)"
  });

  clickArea.on("pointerdown", (event) => {
    if (!metrics.visible) {
      return;
    }

    event.cancelBubble = true;
    const pointer = alignmentStage.getPointerPosition();
    const centeredThumbX = pointer.x - metrics.thumbWidth / 2;
    setAlignmentScrollFromThumbX(centeredThumbX, alignmentLength, sampleCount);
    redrawAlignmentViewer();
  });

  layer.add(clickArea);
}

function renderAlignmentEmpty(message) {
  const subtitle = getAlignmentSubtitle();
  if (subtitle) {
    subtitle.textContent = message;
  }

  alignmentStage.width(getAlignmentStageWidth());
  alignmentStage.height(160);

  alignmentDrawLayer.destroyChildren();
  alignmentInteractionLayer.destroyChildren();

  alignmentDrawLayer.add(new Konva.Rect({
    x: 0,
    y: 0,
    width: getAlignmentStageWidth(),
    height: 160,
    fill: "white",
    listening: false
  }));

  alignmentDrawLayer.add(new Konva.Text({
    x: 14,
    y: 18,
    width: Math.max(1, getAlignmentStageWidth() - 28),
    text: message,
    fontSize: 13,
    fill: "#6b7280",
    listening: false
  }));

  alignmentStage.draw();
  syncSidebarHeightToViewerColumn();
}

function redrawAlignmentViewer() {
  const panel = getAlignmentPanel();
  const blockId = state.activeAlignmentBlockId;
  const alignmentData = getAlignmentData(blockId);
  const containerWidth = getAlignmentStageWidth();

  const stateUnchanged =
    blockId === lastAlignmentRenderState.blockId &&
    state.alignmentFocusedSnpColumn === lastAlignmentRenderState.focusedSnpColumn &&
    state.alignmentZoomX === lastAlignmentRenderState.alignmentZoomX &&
    state.alignmentScrollX === lastAlignmentRenderState.alignmentScrollX &&
    containerWidth === lastAlignmentRenderState.containerWidth;

  if (stateUnchanged) {
    return;
  }

  lastAlignmentRenderState.blockId = blockId;
  lastAlignmentRenderState.focusedSnpColumn = state.alignmentFocusedSnpColumn;
  lastAlignmentRenderState.alignmentZoomX = state.alignmentZoomX;
  lastAlignmentRenderState.alignmentScrollX = state.alignmentScrollX;
  lastAlignmentRenderState.containerWidth = containerWidth;

  if (!blockId) {
    if (panel) {
      panel.classList.add("hidden");
    }
    syncSidebarHeightToViewerColumn();
    return;
  }

  if (panel) {
    panel.classList.remove("hidden");
  }

  if (!alignmentData) {
    renderAlignmentEmpty(`No alignment available for block ${blockId}.`);
    return;
  }

  const sampleNames = getAlignmentSampleNames(alignmentData);
  const alignmentLength = getAlignmentLength(alignmentData);

  if (sampleNames.length === 0 || alignmentLength === 0) {
    renderAlignmentEmpty(`Empty alignment for block ${blockId}.`);
    return;
  }

  normalizeAlignmentScrollX(alignmentLength);
  lastAlignmentRenderState.alignmentScrollX = state.alignmentScrollX;

  const subtitle = getAlignmentSubtitle();
  if (subtitle) {
    subtitle.textContent = `Block ${blockId} (${alignmentLength} bp)`;
  }

  const stageHeight = getAlignmentStageHeight(sampleNames.length);

  alignmentStage.width(getAlignmentStageWidth());
  alignmentStage.height(stageHeight);

  alignmentDrawLayer.destroyChildren();
  alignmentInteractionLayer.destroyChildren();

  alignmentDrawLayer.add(new Konva.Rect({
    x: 0,
    y: 0,
    width: getAlignmentStageWidth(),
    height: stageHeight,
    fill: "white",
    listening: false
  }));

  const snpColumns = getBlockSnpAlignmentColumns(blockId);
  drawAlignmentAxis(alignmentDrawLayer, alignmentLength, snpColumns);
  drawAlignmentRows(
    alignmentDrawLayer,
    alignmentData,
    sampleNames,
    alignmentLength,
    snpColumns
  );
  drawAlignmentScrollbar(alignmentInteractionLayer, alignmentLength, sampleNames.length);

  alignmentStage.draw();
  syncSidebarHeightToViewerColumn();
}

function updateActiveAlignmentViewer() {
  const target = getActiveAlignmentTarget();

  if (!target || !target.blockId) {
    state.activeAlignmentBlockId = null;
    state.alignmentFocusedSnpColumn = null;
    redrawAlignmentViewer();
    return;
  }

  const nextBlockId = target.blockId;
  const nextFocusColumn = target.focusColumn;

  if (state.activeAlignmentBlockId !== nextBlockId) {
    state.activeAlignmentBlockId = nextBlockId;
    state.alignmentScrollX = 0;
    state.alignmentZoomX = 1;
  }

  state.alignmentFocusedSnpColumn = nextFocusColumn;

  const alignmentData = getAlignmentData(nextBlockId);
  const alignmentLength = getAlignmentLength(alignmentData);

  if (nextFocusColumn !== null) {
    centerAlignmentOnColumn(nextFocusColumn, alignmentLength);
  }

  redrawAlignmentViewer();
}

function moveAlignmentByViewportFraction(direction) {
  const alignmentData = getAlignmentData(state.activeAlignmentBlockId);
  const alignmentLength = getAlignmentLength(alignmentData);

  if (alignmentLength === 0) {
    return;
  }

  const stepPx = getAlignmentViewportWidth() * ALIGNMENT.panFraction;
  state.alignmentScrollX = clampAlignmentScrollX(
    state.alignmentScrollX + direction * stepPx,
    alignmentLength
  );

  requestAlignmentRedraw();
}

function zoomAlignmentAroundCenter(nextZoom) {
  const alignmentData = getAlignmentData(state.activeAlignmentBlockId);
  const alignmentLength = getAlignmentLength(alignmentData);

  if (alignmentLength === 0) {
    return;
  }

  const previousCharWidth = getAlignmentCharWidth();
  const centerColumn = (state.alignmentScrollX + getAlignmentViewportWidth() / 2)
    / previousCharWidth;

  state.alignmentZoomX = Math.max(0.5, Math.min(4, nextZoom));

  const nextCharWidth = getAlignmentCharWidth();
  state.alignmentScrollX = centerColumn * nextCharWidth - getAlignmentViewportWidth() / 2;
  normalizeAlignmentScrollX(alignmentLength);
  redrawAlignmentViewer();
}

function startAlignmentViewportDrag(pointerX) {
  const alignmentData = getAlignmentData(state.activeAlignmentBlockId);
  const alignmentLength = getAlignmentLength(alignmentData);

  if (getAlignmentMaxScrollX(alignmentLength) <= 0) {
    return;
  }

  state.isDraggingAlignmentViewport = true;
  state.alignmentDragStartPointerX = pointerX;
  state.alignmentDragStartScrollX = state.alignmentScrollX;
  setBodyCursor("grabbing");
}

function updateAlignmentViewportDrag(pointerX) {
  const alignmentData = getAlignmentData(state.activeAlignmentBlockId);
  const alignmentLength = getAlignmentLength(alignmentData);
  const deltaX = pointerX - state.alignmentDragStartPointerX;

  state.alignmentScrollX = clampAlignmentScrollX(
    state.alignmentDragStartScrollX - deltaX,
    alignmentLength
  );
  requestAlignmentRedraw();
}

function updateAlignmentScrollbarDrag(pointerX) {
  const alignmentData = getAlignmentData(state.activeAlignmentBlockId);
  const alignmentLength = getAlignmentLength(alignmentData);
  const sampleCount = getAlignmentSampleNames(alignmentData).length;

  setAlignmentScrollFromThumbX(
    pointerX - state.alignmentScrollbarDragOffsetX,
    alignmentLength,
    sampleCount
  );
  requestAlignmentRedraw();
}

function stopAlignmentDrag() {
  const wasDragging = state.isDraggingAlignmentViewport
    || state.isDraggingAlignmentScrollbar;

  state.isDraggingAlignmentViewport = false;
  state.isDraggingAlignmentScrollbar = false;
  state.alignmentScrollbarDragOffsetX = 0;

  if (wasDragging) {
    setBodyCursor("default");
    setAlignmentCursor("");
  }
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
  state.isHoveringInteractiveFeature = false;
  _lastResolvedHoverKey = null;

  if (wasDragging) {
    state.suppressHover = false;
    setBodyCursor("default");
    setViewerCursor("");
  }
}

function setupColumnResizer() {
  const contentRow = document.getElementById("content-row");
  const rightColumn = document.getElementById("right-column");
  const resizer = document.getElementById("column-resizer");

  if (!contentRow || !rightColumn || !resizer) {
    return;
  }

  let isResizing = false;

  resizer.addEventListener("pointerdown", (event) => {
    isResizing = true;
    resizer.classList.add("is-dragging");
    setBodyCursor("col-resize");
    resizer.setPointerCapture(event.pointerId);
    showViewerBusyOverlay(" ");
    event.preventDefault();
  });

  resizer.addEventListener("pointermove", (event) => {
    if (!isResizing) {
      return;
    }

    const rowRect = contentRow.getBoundingClientRect();
    const maxSidebarWidth = rowRect.width * CONFIG.sidebarMaxWidthRatio;
    const proposedWidth = rowRect.right - event.clientX;
    const sidebarWidth = Math.max(
      CONFIG.sidebarMinWidth,
      Math.min(maxSidebarWidth, proposedWidth)
    );

    rightColumn.style.flexBasis = `${sidebarWidth}px`;
    rightColumn.style.width = `${sidebarWidth}px`;
  });

  function stopColumnResize(event) {
    if (!isResizing) {
      return;
    }

    isResizing = false;
    resizer.classList.remove("is-dragging");
    setBodyCursor("default");

    if (resizer.hasPointerCapture(event.pointerId)) {
      resizer.releasePointerCapture(event.pointerId);
    }

    normalizeScrollX();
    showViewerBusyOverlay("Rendering viewer\u2026");
    requestStageRedraw();
    requestAlignmentRedraw();
    requestDotplotRedraw();
    syncSidebarHeightToViewerColumn();
    requestAnimationFrame(() => {
      requestAnimationFrame(() => {
        hideViewerBusyOverlay();
      });
    });
  }

  resizer.addEventListener("pointerup", stopColumnResize);
  resizer.addEventListener("pointercancel", stopColumnResize);
}

function syncSidebarHeightToViewerColumn() {
  const viewerColumn = document.getElementById("viewer-column");
  const rightColumn = document.getElementById("right-column");

  if (!viewerColumn || !rightColumn) {
    return;
  }

  const viewerHeight = viewerColumn.getBoundingClientRect().height;
  rightColumn.style.height = `${viewerHeight}px`;
}

function setupFloatingTooltips() {
  const tooltip = document.getElementById("floating-tooltip");

  if (!tooltip) {
    return;
  }

  document.addEventListener("mouseover", (event) => {
    const trigger = event.target.closest(".info-tooltip");

    if (!trigger) {
      return;
    }

    tooltip.textContent = trigger.dataset.tooltip || "";
    tooltip.style.display = "block";

    const triggerRect = trigger.getBoundingClientRect();
    const tooltipRect = tooltip.getBoundingClientRect();
    const margin = 12;

    let left = triggerRect.left + triggerRect.width / 2 - tooltipRect.width / 2;
    let top = triggerRect.top - tooltipRect.height - 8;

    left = Math.max(margin, Math.min(left, window.innerWidth - tooltipRect.width - margin));

    if (top < margin) {
      top = triggerRect.bottom + 8;
    }

    tooltip.style.left = `${left}px`;
    tooltip.style.top = `${top}px`;
  });

  document.addEventListener("mouseout", (event) => {
    const trigger = event.target.closest(".info-tooltip");

    if (!trigger) {
      return;
    }

    tooltip.style.display = "none";
  });
}

setupFloatingTooltips();

function getViewerToolbarHeight() {
  // The toolbar is no longer overlaid on the canvas; it lives above it in normal
  // document flow, so no space needs to be reserved at the top of the Konva stage.
  return 0;
}

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
  state.isHoveringInteractiveFeature = false;
  stopDrag();
  setViewerCursor("");
});

alignmentStage.on("pointerdown", (event) => {
  if (event.target !== alignmentStage) {
    return;
  }

  const pointer = alignmentStage.getPointerPosition();
  if (!pointer) {
    return;
  }

  const alignmentData = getAlignmentData(state.activeAlignmentBlockId);
  const sampleCount = getAlignmentSampleNames(alignmentData).length;
  const scrollbarY = getAlignmentScrollbarY(sampleCount);

  if (pointer.y >= scrollbarY) {
    return;
  }

  startAlignmentViewportDrag(pointer.x);
});

alignmentStage.on("pointermove", () => {
  const pointer = alignmentStage.getPointerPosition();
  if (!pointer) {
    return;
  }

  if (state.isDraggingAlignmentViewport) {
    updateAlignmentViewportDrag(pointer.x);
    return;
  }

  if (state.isDraggingAlignmentScrollbar) {
    updateAlignmentScrollbarDrag(pointer.x);
    return;
  }

  if (state.activeAlignmentBlockId && getAlignmentData(state.activeAlignmentBlockId)) {
    setAlignmentCursor("grab");
  } else {
    setAlignmentCursor("");
  }
});

alignmentStage.on("pointerup", stopAlignmentDrag);
alignmentStage.on("pointerleave", () => {
  stopAlignmentDrag();
  setAlignmentCursor("");
});

function buildSearchIndexes() {
  const {
    blockIdToFeatureId,
    snpKeyToFeatureId,
    featureIdToFeatureType,
    featureIdToRegionRange,
    sampleByName
  } = searchIndexes;

  for (const sample of REGION_DATA.samples) {
    sampleByName.set(sample.sample, sample);

    for (const block of sample.blocks) {
      const numericBlockId = Number(block.block_id);
      const featureId = block.feature_id;

      if (!blockIdToFeatureId.has(numericBlockId)) {
        blockIdToFeatureId.set(numericBlockId, featureId);
      }

      featureIdToFeatureType.set(featureId, "block");

      if (!featureIdToRegionRange.has(featureId)) {
        featureIdToRegionRange.set(featureId, {
          start: block.block_start_in_region,
          end: block.block_end_in_region
        });
      }
    }

    for (const snp of sample.snps) {
      const snpKey = `${snp.block_id}:${snp.aln_pos}`;
      const featureId = snp.feature_id;

      if (!snpKeyToFeatureId.has(snpKey)) {
        snpKeyToFeatureId.set(snpKey, featureId);
      }

      featureIdToFeatureType.set(featureId, "snp");

      if (!featureIdToRegionRange.has(featureId)) {
        featureIdToRegionRange.set(featureId, {
          start: snp.pos_in_region,
          end: snp.pos_in_region
        });
      }
    }
  }
}

function setSearchMode(mode) {
  _searchState.mode = mode;

  const chips = document.querySelectorAll(".search-mode-chip");
  for (const chip of chips) {
    chip.classList.toggle("active", chip.dataset.mode === mode);
  }

  const inputEl = document.getElementById("search-input");
  const sampleSelect = document.getElementById("search-sample-select");
  const statusEl = document.getElementById("search-status");

  if (inputEl) {
    if (mode === "id") {
      inputEl.placeholder = "Block or SNP ID";
    } else if (mode === "region") {
      inputEl.placeholder = "Region coordinate (bp)";
    } else {
      inputEl.placeholder = "Source coordinate (bp)";
    }
  }

  if (sampleSelect) {
    sampleSelect.classList.toggle("hidden", mode !== "source");
  }

  if (statusEl) {
    statusEl.textContent = "";
  }
}

function resolveIdSearch(query) {
  const trimmed = query.trim();

  if (!trimmed) {
    return { error: "Please enter an ID." };
  }

  if (trimmed.includes("::")) {
    const featureType = searchIndexes.featureIdToFeatureType.get(trimmed);
    if (featureType) {
      return { featureId: trimmed, featureType };
    }
    return { error: `Feature not found: "${trimmed}".` };
  }

  if (trimmed.includes(":")) {
    const featureId = searchIndexes.snpKeyToFeatureId.get(trimmed);
    if (featureId) {
      return { featureId, featureType: "snp" };
    }
    return { error: `SNP not found: "${trimmed}".` };
  }

  const blockId = parseInt(trimmed, 10);
  if (!Number.isNaN(blockId) && String(blockId) === trimmed) {
    const featureId = searchIndexes.blockIdToFeatureId.get(blockId);
    if (featureId) {
      return { featureId, featureType: "block" };
    }
    return { error: `Block ${blockId} not found.` };
  }

  return { error: `Unrecognized ID: "${trimmed}".` };
}

function resolveRegionPositionSearch(position) {
  if (!Number.isFinite(position) || position <= 0) {
    return { error: "Please enter a valid position." };
  }

  if (position > REGION_DATA.max_region_length) {
    return { error: `Position out of range (max: ${REGION_DATA.max_region_length}).` };
  }

  return { position };
}

function resolveSourcePositionSearch(sampleName, position) {
  if (!Number.isFinite(position) || position <= 0) {
    return { error: "Please enter a valid position." };
  }

  const sample = searchIndexes.sampleByName.get(sampleName);

  if (!sample) {
    return { error: `Sample "${sampleName}" not found.` };
  }

  const posInRegion = position - sample.region_start_in_source_seq + 1;

  if (posInRegion < 1 || posInRegion > sample.region_length) {
    return {
      error: `Position ${position} is outside the region for "${sampleName}" ` +
        `(region: ${sample.region_start_in_source_seq}–` +
        `${sample.region_start_in_source_seq + sample.region_length - 1}).`
    };
  }

  return { posInRegion };
}

function centerRegionOnRange(start, end) {
  const rangeLength = Math.max(1, end - start + 1);
  const targetVisibleBp = Math.max(CONFIG.targetVisibleBp, rangeLength * 3);
  const desiredZoom = REGION_DATA.max_region_length / targetVisibleBp;
  const newZoom = Math.min(
    getMaxZoomX(),
    Math.max(getInitialZoomX(), desiredZoom)
  );
  state.zoomX = newZoom;
  const centerBp = (start + end) / 2;
  const visibleSpan = getVisibleBpSpan();
  setVisibleStartBp(centerBp - visibleSpan / 2);
  requestStageRedraw();
}

function centerRegionOnPosition(position) {
  const visibleSpan = getVisibleBpSpan();
  setVisibleStartBp(position - visibleSpan / 2);
  requestStageRedraw();
}

function centerRegionOnPositionWithZoom(position) {
  const desiredZoom = REGION_DATA.max_region_length / CONFIG.targetVisibleBp;
  state.zoomX = Math.min(getMaxZoomX(), Math.max(getInitialZoomX(), desiredZoom));
  const visibleSpan = getVisibleBpSpan();
  setVisibleStartBp(position - visibleSpan / 2);
  requestStageRedraw();
}

function runFeatureSearch() {
  const statusEl = document.getElementById("search-status");
  const inputEl = document.getElementById("search-input");
  const sampleSelect = document.getElementById("search-sample-select");

  if (!statusEl || !inputEl) {
    return;
  }

  statusEl.textContent = "";
  const mode = _searchState.mode;
  const rawInput = inputEl.value;

  if (mode === "id") {
    const result = resolveIdSearch(rawInput);
    if (result.error) {
      statusEl.textContent = result.error;
      return;
    }
    const range = searchIndexes.featureIdToRegionRange.get(result.featureId);
    if (range) {
      centerRegionOnRange(range.start, range.end);
    }
    setPinnedFeature(result.featureType, result.featureId);
    return;
  }

  if (mode === "region") {
    const position = Number(rawInput.trim());
    const result = resolveRegionPositionSearch(position);
    if (result.error) {
      statusEl.textContent = result.error;
      return;
    }
    centerRegionOnPositionWithZoom(result.position);
    return;
  }

  if (mode === "source") {
    const sampleName = sampleSelect ? sampleSelect.value : "";
    if (!sampleName) {
      statusEl.textContent = "Please select a sample.";
      return;
    }
    const position = Number(rawInput.trim());
    const result = resolveSourcePositionSearch(sampleName, position);
    if (result.error) {
      statusEl.textContent = result.error;
      return;
    }
    centerRegionOnPositionWithZoom(result.posInRegion);
  }
}

function setupSearchUI() {
  const sampleSelect = document.getElementById("search-sample-select");
  if (sampleSelect) {
    const placeholder = document.createElement("option");
    placeholder.value = "";
    placeholder.textContent = "\u2014 sample \u2014";
    placeholder.disabled = true;
    placeholder.selected = true;
    sampleSelect.appendChild(placeholder);

    for (const sampleName of getSampleOrder()) {
      const option = document.createElement("option");
      option.value = sampleName;
      option.textContent = sampleName;
      sampleSelect.appendChild(option);
    }
  }

  const searchToggle = document.getElementById("search-toggle");
  const searchRow = document.getElementById("search-row");
  if (searchToggle && searchRow) {
    searchToggle.addEventListener("click", () => {
      _searchState.isOpen = !_searchState.isOpen;
      searchRow.classList.toggle("hidden", !_searchState.isOpen);
      searchToggle.classList.toggle("search-open", _searchState.isOpen);
      if (_searchState.isOpen) {
        const input = document.getElementById("search-input");
        if (input) {
          input.focus();
        }
      }
      requestStageRedraw();
    });
  }

  const chips = document.querySelectorAll(".search-mode-chip");
  for (const chip of chips) {
    chip.addEventListener("click", () => {
      setSearchMode(chip.dataset.mode);
    });
  }

  const goButton = document.getElementById("search-go");
  if (goButton) {
    goButton.addEventListener("click", runFeatureSearch);
  }

  const inputEl = document.getElementById("search-input");
  if (inputEl) {
    inputEl.addEventListener("keydown", (event) => {
      if (event.key === "Enter") {
        event.preventDefault();
        runFeatureSearch();
      }
    });
  }

  setSearchMode("id");
}

function initDerivedData() {
  derivedData.sampleOrder = REGION_DATA.samples.map(function(s) {
    return s.sample;
  });

  const trackNameSet = new Set();
  for (const sample of REGION_DATA.samples) {
    for (const track of (sample.gff_tracks || [])) {
      trackNameSet.add(track.track_name);
    }
  }
  derivedData.allGffTrackNames = [...trackNameSet].sort();
  derivedData.gffTrackColorByName = new Map(
    derivedData.allGffTrackNames.map(function(name, index) {
      return [name, GFF_TRACK.colors[index % GFF_TRACK.colors.length]];
    })
  );

  const blockPositions = new Map();
  for (const sample of REGION_DATA.samples) {
    for (const block of sample.blocks) {
      if (!blockPositions.has(block.feature_id)) {
        blockPositions.set(block.feature_id, Number(block.block_start_in_region));
      }
    }
  }
  derivedData.orderedBlockFeatureIds = [...blockPositions.entries()]
    .sort(function(a, b) { return a[1] - b[1]; })
    .map(function(entry) { return entry[0]; });

  const snpPositions = new Map();
  for (const sample of REGION_DATA.samples) {
    for (const snp of sample.snps) {
      if (!snpPositions.has(snp.feature_id)) {
        snpPositions.set(snp.feature_id, Number(snp.pos_in_region));
      }
    }
  }
  derivedData.orderedSnpFeatureIds = [...snpPositions.entries()]
    .sort(function(a, b) { return a[1] - b[1]; })
    .map(function(entry) { return entry[0]; });

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

function initializeViewer() {
  state.featureGroups = buildFeatureGroups(REGION_DATA);
  initDerivedData();
  buildSearchIndexes();
  renderAnalysisSettings();
  renderSidebarDefault();

  state.zoomX = getInitialZoomX();
  state.scrollX = 0;

  setupColumnResizer();
  setupSearchUI();
  setupModeSwitch();
  updateViewerModeInfoTooltip("browser");
  buildDotplotPairIndexes();
  setupDotplotUI();
  setupWheelScrolling();

  redrawStage();
  redrawAlignmentViewer();
  syncSidebarHeightToViewerColumn();

  requestAnimationFrame(() => {
    hideRenderingOverlay();
  });
}

showRenderingOverlay();

requestAnimationFrame(() => {
  initializeViewer();
});

window.addEventListener("resize", () => {
  normalizeScrollX();
  requestStageRedraw();
  requestAlignmentRedraw();
  syncSidebarHeightToViewerColumn();
  requestDotplotRedraw();
});

document.getElementById("feature-prev").addEventListener("click", () => {
  pinNeighborFeature(-1);
});

document.getElementById("feature-next").addEventListener("click", () => {
  pinNeighborFeature(1);
});

document.getElementById("feature-center").addEventListener("click", () => {
  centerPinnedFeature();
});

document.getElementById("zoom-in").addEventListener("click", () => {
  if (isDotplotModeActive()) {
    _dotplotState.zoom = Math.min(DOTPLOT_ZOOM_MAX, _dotplotState.zoom * DOTPLOT_ZOOM_STEP);
    requestDotplotRedraw();
    return;
  }

  zoomAroundViewportCenter(state.zoomX * getZoomFactor());
  redrawStage();
});

document.getElementById("zoom-out").addEventListener("click", () => {
  if (isDotplotModeActive()) {
    _dotplotState.zoom = Math.max(DOTPLOT_ZOOM_MIN, _dotplotState.zoom / DOTPLOT_ZOOM_STEP);
    requestDotplotRedraw();
    return;
  }

  zoomAroundViewportCenter(state.zoomX / getZoomFactor());
  redrawStage();
});

document.getElementById("zoom-reset").addEventListener("click", () => {
  resetActiveViewerZoom();
});


function getWheelDeltaX(event) {
  if (Math.abs(event.deltaX) > Math.abs(event.deltaY)) {
    return event.deltaX;
  }

  if (event.shiftKey) {
    return event.deltaY;
  }

  return 0;
}

function setActiveKeyboardViewer(viewerName) {
  state.activeKeyboardViewer = viewerName;
}

function setupWheelScrolling() {
  const viewerElement = document.getElementById("viewer");
  const alignmentElement = document.getElementById("alignment-viewer");

  viewerElement.addEventListener("wheel", (event) => {
    const deltaX = getWheelDeltaX(event);

    if (deltaX === 0) {
      return;
    }

    event.preventDefault();
    state.scrollX = clampScrollX(state.scrollX + deltaX);
    requestStageRedraw();
  }, { passive: false });

  alignmentElement.addEventListener("wheel", (event) => {
    const deltaX = getWheelDeltaX(event);

    if (deltaX === 0) {
      return;
    }

    event.preventDefault();
    const alignmentData = getAlignmentData(state.activeAlignmentBlockId);
    const alignmentLength = getAlignmentLength(alignmentData);
    state.alignmentScrollX = clampAlignmentScrollX(
      state.alignmentScrollX + deltaX,
      alignmentLength
    );
    requestAlignmentRedraw();
  }, { passive: false });

  viewerElement.addEventListener("pointerdown", () => {
    setActiveKeyboardViewer("region");
  });

  alignmentElement.addEventListener("pointerdown", () => {
    setActiveKeyboardViewer("alignment");
  });
}

window.addEventListener("keydown", (event) => {
  const activeElement = document.activeElement;
  const isSearchInputFocused = activeElement && activeElement.id === "search-input";
  const key = event.key.toLowerCase();

  if (isSearchInputFocused) {
    return;
  }

  if (key === "r" && !event.ctrlKey && !event.metaKey) {
    event.preventDefault();
    resetActiveViewerZoom();
    return;
  }

  if (key === "f") {
    event.preventDefault();
    centerPinnedFeature();
    return;
  }

  if (event.key === "Tab") {
    const hasPinnedFeature = Boolean(state.pinnedFeatureId && state.pinnedFeatureType);

    if (!hasPinnedFeature) {
      return;
    }

    event.preventDefault();

    if (event.shiftKey) {
      pinNeighborFeature(-1);
    } else {
      pinNeighborFeature(1);
    }

    return;
  }

  if (event.key !== "ArrowLeft" && event.key !== "ArrowRight") {
    return;
  }

  event.preventDefault();
  const direction = event.key === "ArrowLeft" ? -1 : 1;

  if (isDotplotModeActive()) {
    moveDotplotByViewportFraction(direction, 0.1);
    return;
  }

  if (state.activeKeyboardViewer === "alignment") {
    moveAlignmentByViewportFraction(direction);
    return;
  }

  moveByViewportFraction(direction, 0.1);
});

document.getElementById("alignment-snp-prev").addEventListener("click", () => {
  focusNeighborAlignmentSnp(-1);
});

document.getElementById("alignment-snp-next").addEventListener("click", () => {
  focusNeighborAlignmentSnp(1);
});

document.getElementById("alignment-zoom-in").addEventListener("click", () => {
  zoomAlignmentAroundCenter(state.alignmentZoomX * 1.25);
});

document.getElementById("alignment-zoom-out").addEventListener("click", () => {
  zoomAlignmentAroundCenter(state.alignmentZoomX / 1.25);
});

document.getElementById("alignment-zoom-reset").addEventListener("click", () => {
  state.alignmentZoomX = 1;
  state.alignmentScrollX = 0;
  redrawAlignmentViewer();
});
