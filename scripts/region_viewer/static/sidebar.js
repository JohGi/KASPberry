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
