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
