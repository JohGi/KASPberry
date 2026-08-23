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
    return `Kimura 2-parameter distances are computed from unmasked MAFFT block alignments with EMBOSS \`distmat\` using \`-nucmethod 2\`.\n\nEMBOSS reports distances per 100 sites. KASPberry converts these values to substitutions per site.\n\nReferences: Kimura 1980; Rice, Longden & Bleasby 2000.`;
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
  browser: {
    snps: [
      "This view shows locally collinear blocks shared by all input genotypes and SNPs detected within these blocks. Diagnostic SNP candidates are shown in red. Rejected candidates are shown in gray.",
      "",
      "Locally collinear blocks are identified with SibeliaZ. Blocks are retained if they occur once in each genotype, have the same strand orientation, and pass the configured minimum block length.",
      "",
      "Simple repeats and low-complexity regions are hard-masked with RepeatMasker. If a custom repeat library is provided, matching sequences are also masked. Masked block sequences are aligned with MAFFT. Biallelic SNPs are detected with SeqTUI using the configured minimum SNP flank, i.e. the minimum number of perfectly aligned bases required on both sides of each SNP, with no indel, additional SNP, or N.",
      "",
      "Diagnostic SNPs are retained if all genotypes in each group have the same allele and the two groups have different alleles.",
      "",
      "References: SibeliaZ, Minkin & Medvedev 2020; RepeatMasker, Tarailo-Graovac & Chen 2009; MAFFT, Katoh & Standley 2013; SeqTUI, Ranwez 2026."
    ].join("\n"),
    kasp: [
      "This view shows locally collinear blocks shared by all input genotypes and SNPs detected within these blocks. Final KASP SNP candidates are shown in red. Rejected candidates are shown in gray.",
      "",
      "Locally collinear blocks are identified with SibeliaZ. Blocks are retained if they occur once in each genotype, have the same strand orientation, and pass the configured minimum block length.",
      "",
      "Simple repeats and low-complexity regions are hard-masked with RepeatMasker. If a custom repeat library is provided, matching sequences are also masked. Masked block sequences are aligned with MAFFT. Biallelic SNPs are detected with SeqTUI using the configured minimum SNP flank, i.e. the minimum number of perfectly aligned bases required on both sides of each SNP, with no indel, additional SNP, or N.",
      "",
      "Diagnostic SNPs are selected from the configured genotype groups. KASP assays are designed with PolyMarker. Assays are tested with MFEprimer for in-silico amplification specificity, dimers, and hairpins. SNPs with at least one assay that passes in-silico screening are retained as final candidates.",
      "",
      "References: SibeliaZ, Minkin & Medvedev 2020; RepeatMasker, Tarailo-Graovac & Chen 2009; MAFFT, Katoh & Standley 2013; SeqTUI, Ranwez 2026; PolyMarker, Ramirez-Gonzalez et al. 2015; MFEprimer, Wang et al. 2019."
    ].join("\n")
  },

  dotplot: [
    "This view shows pairwise dotplots between input regional sequences.",
    "",
    "For each selected genotype pair, sequences are aligned with Minimap2 using the `asm5` preset. The alignment is displayed as a dotplot with blastn2dotplots.",
    "",
    "References: Minimap2, Li 2018; blastn2dotplots, Okuno et al. 2025."
  ].join("\n")
};


function updateViewerModeInfoTooltip(mode) {
  const infoTooltip = document.getElementById("viewer-mode-info-tooltip");

  if (!infoTooltip) {
    return;
  }

  infoTooltip.dataset.tooltip = mode === "browser"
    ? VIEWER_MODE_INFO_TOOLTIPS.browser[REGION_DATA.mode] || ""
    : VIEWER_MODE_INFO_TOOLTIPS[mode] || "";
}

function renderDistanceMatrix(matrix, { embedded = false, showTitle = true } = {}) {
  const containerClass = embedded ? "distance-matrix-content" : "distance-matrix-card";

  if (!matrix || !matrix.labels || !matrix.values) {
    return `
      <div class="${containerClass}">
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
    <div class="${containerClass}">
      ${showTitle ? `
        <p class="distance-matrix-title">
          ${escapeHtml(matrix.title || "Distance matrix")}
          ${renderInfoTooltip(tooltip)}
        </p>
      ` : ""}
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

function formatSummaryCount(value) {
  const numericValue = Number(value);
  return Number.isFinite(numericValue) ? formatNumber(Math.round(numericValue)) : "NA";
}

function formatSummaryBp(value) {
  const numericValue = Number(value);
  return Number.isFinite(numericValue) ? `${formatNumber(numericValue, 2)} bp` : "NA";
}

function formatSummaryPercent(value) {
  const numericValue = Number(value);
  return Number.isFinite(numericValue) ? `${formatNumber(numericValue, 2)}%` : "NA";
}

function renderSummaryRows(entries) {
  let html = `<div class="kv">`;

  for (const [label, value] of entries) {
    html += `<div class="key">${escapeHtml(label)}</div><div>${escapeHtml(String(value))}</div>`;
  }

  return `${html}</div>`;
}

function renderSummarySection(title, entries, emphasized = false) {
  return `
    <section class="analysis-summary-section${emphasized ? " analysis-summary-section-emphasized" : ""}">
      <h4>${escapeHtml(title)}</h4>
      ${renderSummaryRows(entries)}
    </section>
  `;
}

function subtractSummaryCounts(total, retained) {
  const totalCount = Number(total);
  const retainedCount = Number(retained);

  if (!Number.isFinite(totalCount) || !Number.isFinite(retainedCount)) {
    return "NA";
  }

  return formatSummaryCount(totalCount - retainedCount);
}

function renderTreeResult(label, value, children = "") {
  return `
    <li class="analysis-summary-tree-result${children ? " analysis-summary-tree-result-with-children" : ""}">
      <div class="analysis-summary-tree-row">
        <span class="analysis-summary-tree-label">${escapeHtml(label)}</span>
        <span class="analysis-summary-tree-value">${escapeHtml(String(value))}</span>
      </div>
      ${children ? `
        <ul class="analysis-summary-tree-children">
          ${children}
        </ul>
      ` : ""}
    </li>
  `;
}

function renderTreeStage(label, children) {
  return `
    <li class="analysis-summary-tree-stage">
      <span class="analysis-summary-tree-stage-label">${escapeHtml(label)}</span>
      <ul class="analysis-summary-tree-children">
        ${children}
      </ul>
    </li>
  `;
}

const ASSAY_FAILURE_REASON_ENTRIES = [
  ["missing_target_amplicon", "Missing target amplicon"],
  ["multiple_expected_amplicons", "Multiple expected amplicons"],
  ["unexpected_allele_amplicon", "Unexpected allele amplicon"],
  ["noncanonical_amplicon", "Noncanonical amplicon"],
  ["dimer", "Dimer"],
  ["hairpin", "Hairpin"]
];

function renderAssayFailureReasonSummary(reasons) {
  return `
    <details class="analysis-summary-reasons">
      <summary>Failure reasons</summary>
      ${renderSummaryRows(ASSAY_FAILURE_REASON_ENTRIES.map(([reason, label]) => [
        label,
        formatSummaryCount(reasons?.[reason] ?? 0)
      ]))}
    </details>
  `;
}

function renderSnpSummary(snpDiscovery, assayDesign, finalCandidates) {
  if (!snpDiscovery) {
    return "";
  }

  const diagnosticSnps = snpDiscovery.diagnostic_snps;
  let diagnosticSnpResult = renderTreeResult(
    "Diagnostic SNPs",
    formatSummaryCount(diagnosticSnps)
  );

  if (REGION_DATA.mode === "kasp" && assayDesign && finalCandidates) {
    const snpsWithAssay = assayDesign.snps_with_assay_proposed;
    const snpsWithPassingAssay = finalCandidates.candidate_snps;
    const assayDesignTree = renderTreeStage("KASP assay design", [
      renderTreeResult(
        "SNPs with no assay proposed",
        subtractSummaryCounts(diagnosticSnps, snpsWithAssay)
      ),
      renderTreeResult(
        "SNPs with assay proposed",
        formatSummaryCount(snpsWithAssay),
        renderTreeStage("In-silico screening", [
          renderTreeResult(
            "SNPs with no passing assay",
            subtractSummaryCounts(snpsWithAssay, snpsWithPassingAssay)
          ),
          renderTreeResult(
            "SNPs with \u22651 passing assay",
            formatSummaryCount(snpsWithPassingAssay)
          )
        ].join(""))
      )
    ].join(""));

    diagnosticSnpResult = renderTreeResult(
      "Diagnostic SNPs",
      formatSummaryCount(diagnosticSnps),
      assayDesignTree
    );
  }

  const diagnosticFilteringTree = renderTreeStage("Diagnostic filtering", [
      renderTreeResult(
        "Non-diagnostic SNPs",
        formatSummaryCount(snpDiscovery.non_diagnostic_snps)
      ),
      diagnosticSnpResult
    ].join(""));

  return `
    <section class="analysis-summary-section analysis-summary-tree-section">
      <h4>SNPs</h4>
      <ul class="analysis-summary-tree">
        ${renderTreeResult(
          "Detected SNPs",
          formatSummaryCount(snpDiscovery.detected_snps),
          diagnosticFilteringTree
        )}
      </ul>
    </section>
  `;
}

function renderAssaySummary(assayDesign, validation, failureReasons) {
  if (!assayDesign || !validation) {
    return "";
  }

  return `
    <section class="analysis-summary-section analysis-summary-tree-section">
      <h4>KASP assays</h4>
      <ul class="analysis-summary-tree">
        ${renderTreeResult(
          "Proposed assays",
          formatSummaryCount(assayDesign.assays_proposed),
          renderTreeStage("In-silico screening", [
            renderTreeResult(
              "Passing assays",
              formatSummaryCount(validation.assays_passing_validation)
            ),
            renderTreeResult(
              "Failing assays",
              formatSummaryCount(validation.assays_failing_validation)
            )
          ].join(""))
        )}
      </ul>
      ${renderAssayFailureReasonSummary(failureReasons)}
    </section>
  `;
}

function nonEmptyReasonEntries(reasons) {
  return Object.entries(reasons || {}).filter(([, count]) => Number(count) > 0);
}

function renderReasonSummary(title, reasons) {
  const entries = nonEmptyReasonEntries(reasons);

  if (entries.length === 0) {
    return "";
  }

  return `
    <details class="analysis-summary-reasons">
      <summary>${escapeHtml(title)}</summary>
      ${renderSummaryRows(entries.map(([reason, count]) => [
        humanizeSnpWorkflowValue(reason),
        formatSummaryCount(count)
      ]))}
    </details>
  `;
}

function renderMashSummary() {
  return `
    <details class="analysis-summary-mash">
      <summary>Mash distances, whole region</summary>
      ${renderDistanceMatrix(REGION_DATA.mash_matrix, { embedded: true, showTitle: false })}
    </details>
  `;
}

function renderAnalysisSummary() {
  const summaryStats = REGION_DATA.summary_stats || {};
  const sections = [];
  const input = summaryStats.input;
  const globalStats = summaryStats.global;
  const snpDiscovery = summaryStats.snp_discovery;

  if (input && Object.prototype.hasOwnProperty.call(input, "n_genotypes")) {
    sections.push(renderSummarySection("Input", [
      ["Genotypes", formatSummaryCount(input.n_genotypes)]
    ]));
  }

  if (globalStats) {
    sections.push(renderSummarySection("Collinear blocks", [
      ["Retained blocks", formatSummaryCount(globalStats.n_blocks_kept)],
      ["Shortest block", formatSummaryBp(globalStats.min_block_len_bp)],
      ["Longest block", formatSummaryBp(globalStats.max_block_len_bp)],
      ["Mean block length", formatSummaryBp(globalStats.mean_block_len_bp)]
    ]));
  }

  if (snpDiscovery) {
    sections.push(renderSnpSummary(
      snpDiscovery,
      summaryStats.kasp_assay_design,
      summaryStats.final_candidates
    ));
  }

  if (REGION_DATA.mode === "kasp") {
    const assayDesign = summaryStats.kasp_assay_design;
    const validation = summaryStats.in_silico_validation;
    const failureReasons = summaryStats.failure_reasons || {};
    sections.push(renderAssaySummary(assayDesign, validation, failureReasons.assays));
  }

  sections.push(renderMashSummary());

  const content = sections.filter(Boolean).join("");
  return `
    <div class="summary-card analysis-summary-card">
      ${content || '<p class="hint">No summary statistics available.</p>'}
    </div>
  `;
}

function renderAnalysisSettingsRows(entries) {
  let html = `<div class="kv">`;

  for (const [label, value] of entries) {
    html += `<div class="key">${escapeHtml(label)}</div><div>${escapeHtml(String(value))}</div>`;
  }

  return `${html}</div>`;
}

function renderAnalysisSettingsSection(title, content) {
  return `
    <section class="analysis-settings-section">
      <h3>${escapeHtml(title)}</h3>
      ${content}
    </section>
  `;
}

function renderDiagnosticGroups(groups) {
  const entries = Object.entries(groups || {});

  if (entries.length === 0) {
    return renderAnalysisSettingsSection(
      "Diagnostic groups",
      `<p class="hint">No genotype groups available.</p>`
    );
  }

  const groupMarkup = entries.map(([group, genotypes]) => `
    <div class="analysis-settings-group">
      <div class="analysis-settings-group-name">${escapeHtml(group)}</div>
      <div>${escapeHtml((genotypes || []).join(", "))}</div>
    </div>
  `).join("");

  return renderAnalysisSettingsSection(
    "Diagnostic groups",
    `<div class="analysis-settings-groups">${groupMarkup}</div>`
  );
}

function renderAdvancedOptions(settings) {
  const labels = {
    sibeliaz: "SibeliaZ",
    mafft: "MAFFT",
    mfeprimer_specificity: "MFEprimer specificity",
    mfeprimer_dimer: "MFEprimer dimer",
    mfeprimer_hairpin: "MFEprimer hairpin"
  };
  const allowedOptions = REGION_DATA.mode === "kasp"
    ? Object.keys(labels)
    : ["sibeliaz", "mafft"];
  const optionSets = Object.entries(settings.advanced_options || {})
    .filter(([name, options]) => allowedOptions.includes(name) && Array.isArray(options) && options.length > 0);

  if (optionSets.length === 0) {
    return "";
  }

  const optionsMarkup = optionSets.map(([name, options]) => `
    <div class="analysis-settings-option-set">
      <div class="analysis-settings-group-name">${escapeHtml(labels[name])}</div>
      <ul>
        ${options.map(option => `<li><code>${escapeHtml(String(option))}</code></li>`).join("")}
      </ul>
    </div>
  `).join("");

  return renderAnalysisSettingsSection(
    "Advanced options",
    `<div class="analysis-settings-options">${optionsMarkup}</div>`
  );
}

function formatCustomRepeatLibrary(library) {
  if (library === null || library === undefined || String(library).trim() === "") {
    return "None";
  }

  const configuredLibrary = String(library).trim();

  if (/(^|\/)repeatmasker_placeholder\.fa$/.test(configuredLibrary.replaceAll("\\", "/"))) {
    return "None";
  }

  const pathParts = configuredLibrary
    .replaceAll("\\", "/")
    .split("/")
    .filter(Boolean);

  return pathParts[pathParts.length - 1] || configuredLibrary;
}

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

  let html = `<div class="analysis-settings">`;
  html += renderAnalysisSettingsSection("SNP discovery", renderAnalysisSettingsRows([
    ["Minimum block length", minBlock],
    ["Minimum SNP flank", minFlank]
  ]));
  html += renderDiagnosticGroups(settings.snp_groups);

  const repeatMasking = settings.repeat_masking || {};
  html += renderAnalysisSettingsSection("Repeat masking", renderAnalysisSettingsRows([
    ["Simple repeats and low complexity", "Yes"],
    ["Custom repeat library", formatCustomRepeatLibrary(repeatMasking.custom_repeat_library)]
  ]));

  if (REGION_DATA.mode === "kasp") {
    const kasp = settings.kasp_assay_design || {};
    html += renderAnalysisSettingsSection("KASP assay design", renderAnalysisSettingsRows([
      ["PolyMarker subgenomes", kasp.polymarker_subgenomes ?? "NA"],
      ["PolyMarker / specificity genotypes", (kasp.genotypes || []).join(", ") || "NA"],
      ["MFEprimer minimum binding Tm", kasp.mfeprimer_min_tm ?? "NA"],
      ["MFEprimer dimer ΔG cutoff", kasp.mfeprimer_dimer_max_dg ?? "NA"]
    ]));
  }

  html += renderAdvancedOptions(settings);
  html += `</div>`;

  container.innerHTML = html;
}

function renderSampleRegionStats() {
  const sampleStats = REGION_DATA.summary_stats?.samples;

  if (!sampleStats || Object.keys(sampleStats).length === 0) {
    return `
      <div class="summary-card">
        <h3>Per-genotype statistics</h3>
        <p class="hint">No per-sample region statistics available.</p>
      </div>
    `;
  }

  let html = `<div class="summary-card"><h3>Per-genotype statistics</h3>`;

  for (const sampleName of getSampleOrder()) {
    const stats = sampleStats[sampleName];
    html += `<div class="sample-card"><h3>${escapeHtml(sampleName)}</h3>`;

    if (!stats) {
      html += `<p class="hint">No data.</p>`;
    } else {
      const entries = [
        ["Region length", formatSummaryBp(stats.region_length_bp)],
        ["Covered by collinear blocks", formatSummaryPercent(stats.covered_pct_of_region)],
        ["Repeat-masked bases in collinear blocks", formatSummaryPercent(stats.repeat_masked_pct_of_collinear_blocks)]
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
      <h2>Analysis summary</h2>
    </div>
    <div class="sidebar-section">
      ${renderAnalysisSummary()}
    </div>
    <div class="sidebar-section">
      ${renderSampleRegionStats()}
    </div>
  `;
}

function renderSidebarHeader(title, isPinned) {
  return `
    <div class="sidebar-header">
      <div class="sidebar-title-row">
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

const SNP_WORKFLOW_DISPLAY_VALUES = {
  non_diagnostic_allele_pattern: "Non-diagnostic allele pattern",
  no_polymarker_assay: "No PolyMarker assay",
  no_assay_passed_in_silico_validation: "No assay passed in-silico screening"
};

function humanizeSnpWorkflowValue(value) {
  if (value === null || value === undefined || String(value).trim() === "") {
    return "";
  }

  const rawValue = String(value).trim();
  if (SNP_WORKFLOW_DISPLAY_VALUES[rawValue]) {
    return SNP_WORKFLOW_DISPLAY_VALUES[rawValue];
  }

  const normalized = rawValue
    .replace(/[_-]+/g, " ")
    .replace(/\s+/g, " ");

  return normalized.charAt(0).toUpperCase() + normalized.slice(1);
}

function getSnpWorkflowStatus(status, isFinalStage = false) {
  const normalized = String(status || "NOT_RUN").trim().toUpperCase();

  if (normalized === "PASS") {
    return { label: "PASS", className: "pass" };
  }

  if (normalized === "FAIL") {
    return {
      label: isFinalStage ? "REJECTED" : "FAIL",
      className: isFinalStage ? "rejected" : "fail"
    };
  }

  return { label: "Not run", className: "not-run" };
}

function renderSnpWorkflowBadge(status, isFinalStage = false) {
  const display = getSnpWorkflowStatus(status, isFinalStage);
  return `<span class="snp-status-badge ${display.className}">${display.label}</span>`;
}

function renderSnpWorkflowReason(reason) {
  const displayReason = humanizeSnpWorkflowValue(reason);
  if (!displayReason) {
    return "";
  }

  return `<p class="snp-workflow-reason">${escapeHtml(displayReason)}</p>`;
}

function renderSnpWorkflowStage(label, status, reason, isFinalStage = false) {
  const display = getSnpWorkflowStatus(status, isFinalStage);
  const renderedReason = display.className === "not-run"
    ? ""
    : renderSnpWorkflowReason(reason);

  return `
    <div class="snp-workflow-stage">
      <div class="snp-workflow-stage-header">
        <span>${escapeHtml(label)}</span>
        ${renderSnpWorkflowBadge(status, isFinalStage)}
      </div>
      ${renderedReason}
    </div>
  `;
}

function renderSnpWorkflowPipeline(snpResult) {
  const isKaspMode = REGION_DATA.mode === "kasp";
  let stages = renderSnpWorkflowStage(
    "Diagnostic selection",
    snpResult?.diagnostic_status,
    snpResult?.diagnostic_failure_reason
  );

  if (isKaspMode) {
    stages += renderSnpWorkflowStage(
      "PolyMarker design",
      snpResult?.design_status,
      snpResult?.design_failure_reason
    );
    stages += renderSnpWorkflowStage(
      "In-silico screening",
      snpResult?.validation_status,
      snpResult?.validation_failure_reason
    );
  }

  return `
    <div class="sidebar-section snp-workflow-section">
      <h3>Candidate processing</h3>
      <div class="snp-workflow-list">${stages}</div>
    </div>
  `;
}

function renderAssayDetail(label, value, className = "") {
  if (value === null || value === undefined || String(value).trim() === "") {
    return "";
  }

  const valueTag = className === "primer-sequence" ? "code" : "div";
  return `
    <div class="assay-detail-row">
      <div class="key">${escapeHtml(label)}</div>
      <${valueTag} class="${className}">${escapeHtml(String(value))}</${valueTag}>
    </div>
  `;
}

function renderAssayFailureReason(status, reason) {
  if (getSnpWorkflowStatus(status).className !== "fail") {
    return "";
  }

  const displayReason = humanizeSnpWorkflowValue(reason);
  if (!displayReason) {
    return "";
  }

  return `<div class="assay-failure-reason"><span class="key">Failure reason:</span> ${escapeHtml(displayReason)}</div>`;
}

function renderAlleleSpecificPrimer(label, primer, primerWithTail, allele) {
  const normalizedPrimer = String(primer || "").toUpperCase();
  const normalizedWithTail = String(primerWithTail || "").toUpperCase();
  const displayedSequence = normalizedWithTail || normalizedPrimer;

  if (!displayedSequence) {
    return "";
  }

  let renderedSequence = escapeHtml(displayedSequence);
  if (normalizedPrimer && normalizedWithTail.endsWith(normalizedPrimer)) {
    const tail = normalizedWithTail.slice(0, -normalizedPrimer.length);
    const body = normalizedPrimer.slice(0, -1);
    const terminalBase = normalizedPrimer.slice(-1);
    const expectedAllele = String(allele || "").toUpperCase();
    const renderedTerminalBase = terminalBase === expectedAllele
      ? `<strong class="primer-terminal-base">${escapeHtml(terminalBase)}</strong>`
      : escapeHtml(terminalBase);

    renderedSequence = `<span class="primer-tail">${escapeHtml(tail)}</span>${escapeHtml(body)}${renderedTerminalBase}`;
  }

  return `
    <div class="assay-detail-row">
      <div class="key">${escapeHtml(label)}</div>
      <code class="primer-sequence">${renderedSequence}</code>
    </div>
  `;
}

function renderKaspAssays(featureId) {
  if (REGION_DATA.mode !== "kasp") {
    return "";
  }

  const assays = REGION_DATA.assays_by_snp?.[featureId] || [];
  let content = "<p class=\"hint\">No PolyMarker assay available for this SNP.</p>";

  if (assays.length > 0) {
    content = '<div class="snp-assay-list">';

    for (const assay of assays) {
      const assayId = assay.assay_id || "Assay";
      const validationStatus = assay.validation_status;
      content += `
        <details class="snp-assay-card">
          <summary>
            <span>${escapeHtml(assayId)}</span>
            ${renderSnpWorkflowBadge(validationStatus)}
          </summary>
          <div class="snp-assay-details">
            ${renderAssayFailureReason(validationStatus, assay.validation_failure_reason)}
            ${renderAlleleSpecificPrimer(`Allele ${String(assay.first_allele || "").toUpperCase()} primer`, assay.first_primer, assay.first_primer_with_tail, assay.first_allele)}
            ${renderAlleleSpecificPrimer(`Allele ${String(assay.second_allele || "").toUpperCase()} primer`, assay.second_primer, assay.second_primer_with_tail, assay.second_allele)}
            ${renderAssayDetail("Common primer", String(assay.common_primer || "").toUpperCase(), "primer-sequence")}
          </div>
        </details>
      `;
    }

    content += "</div>";
  }

  return `
    <div class="sidebar-section snp-assay-section">
      <h3>KASP assays</h3>
      ${content}
    </div>
  `;
}

function attachSnpSidebarHandlers(blockId) {
  const blockLink = document.getElementById("snp-block-link");

  if (!blockLink) {
    return;
  }

  blockLink.addEventListener("click", () => {
    const featureId = `block::${blockId}`;
    const range = searchIndexes.featureIdToRegionRange.get(featureId);

    if (!range) {
      return;
    }

    state.hoveredFeatureType = null;
    state.hoveredFeatureId = null;
    _lastResolvedHoverKey = null;
    _lastResolvedDotplotHoverKey = null;
    setPinnedFeature("block", featureId);
    centerRegionOnRange(range.start, range.end);
  });
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
  const snpResult = getSnpResult(featureId);

  let html = `
    ${renderSidebarHeader("SNP", isPinned)}
    <p class="hint"><b>ID:</b> ${escapeHtml(title)}</p>
    <p class="hint snp-block-metadata"><b>Block:</b> <button id="snp-block-link" class="snp-block-link" type="button">${escapeHtml(firstInfo.block_id)}</button></p>
    <div class="snp-candidate-summary">
      <span>Candidate status</span>
      ${renderSnpWorkflowBadge(snpResult?.final_status, true)}
      ${renderSnpWorkflowReason(snpResult?.final_failure_reason)}
    </div>
    ${renderSnpWorkflowPipeline(snpResult)}
    ${renderKaspAssays(featureId)}
  `;

  const observationsOpen = REGION_DATA.mode === "kasp" ? "" : " open";
  html += `
    <div class="sidebar-section snp-observations-section">
      <details class="snp-observations"${observationsOpen}>
        <summary>Alleles and positions</summary>
        <div class="snp-observation-cards">
  `;

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

  html += `
        </div>
      </details>
    </div>
  `;

  sidebar.innerHTML = html;
  attachSidebarUnpinHandler();
  attachSnpSidebarHandlers(firstInfo.block_id);
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
