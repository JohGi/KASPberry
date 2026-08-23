function syncSidebarHeightToViewerColumn() {
  const viewerColumn = document.getElementById("viewer-column");
  const rightColumn = document.getElementById("right-column");

  if (!viewerColumn || !rightColumn) {
    return;
  }

  const viewerHeight = viewerColumn.getBoundingClientRect().height;
  rightColumn.style.height = `${viewerHeight}px`;
}

const MAX_SIDEBAR_WIDTH_RATIO = 0.7;

function setupColumnResizer() {
  const contentRow = document.getElementById("content-row");
  const rightColumn = document.getElementById("right-column");
  const resizer = document.getElementById("column-resizer");

  if (!contentRow || !rightColumn || !resizer) {
    return;
  }

  const sidebarMinWidth = Number.parseFloat(
    window.getComputedStyle(rightColumn).minWidth
  );
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
    const maxSidebarWidth = rowRect.width * MAX_SIDEBAR_WIDTH_RATIO;
    const proposedWidth = rowRect.right - event.clientX;
    const sidebarWidth = Math.max(
      sidebarMinWidth,
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
    preserveDotplotViewportForNextRedraw();
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

function setupRejectedSnpToggle() {
  const toggle = document.getElementById("show-rejected-snps");

  if (!toggle) {
    return;
  }

  toggle.checked = state.showRejectedSnps;
  toggle.addEventListener("change", () => {
    state.showRejectedSnps = toggle.checked;
    _hoverIndexDirty = true;
    _lastResolvedHoverKey = null;
    _dotplotHoverIndexDirty = true;
    _lastResolvedDotplotHoverKey = null;

    if (
      state.hoveredFeatureType === "snp" &&
      state.hoveredFeatureId &&
      isRejectedSnp(state.hoveredFeatureId)
    ) {
      clearHoveredFeature();
    }

    if (
      state.pinnedFeatureType === "snp" &&
      state.pinnedFeatureId &&
      isRejectedSnp(state.pinnedFeatureId)
    ) {
      clearPinnedFeature();
    }

    requestStageRedraw();
    if (isDotplotModeActive()) {
      requestDotplotRedraw();
    }
  });
}

setupRejectedSnpToggle();

function getViewerToolbarHeight() {
  // The toolbar is no longer overlaid on the canvas; it lives above it in normal
  // document flow, so no space needs to be reserved at the top of the Konva stage.
  return 0;
}

function setSearchMode(mode) {
  if (isDotplotModeActive() && mode !== "id") {
    mode = "id";
  }

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
    sampleSelect.classList.toggle("hidden", mode !== "source" || isDotplotModeActive());
  }

  if (statusEl) {
    statusEl.textContent = "";
  }
}

function updateSearchModeAvailability(isBrowser) {
  const chips = document.querySelectorAll(".search-mode-chip");
  for (const chip of chips) {
    chip.classList.toggle("hidden", !isBrowser && chip.dataset.mode !== "id");
  }

  if (!isBrowser && _searchState.mode !== "id") {
    setSearchMode("id");
    return;
  }

  const sampleSelect = document.getElementById("search-sample-select");
  if (sampleSelect) {
    sampleSelect.classList.toggle("hidden", !isBrowser || _searchState.mode !== "source");
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
  const targetVisibleBp = Math.max(BROWSER_ZOOM.targetVisibleBp, rangeLength * 3);
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

function centerRegionOnPositionWithZoom(position) {
  const desiredZoom = REGION_DATA.max_region_length / BROWSER_ZOOM.targetVisibleBp;
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
