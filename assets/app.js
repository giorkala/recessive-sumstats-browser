const DATA_ROOT = "data";
const MAX_TABLE_ROWS = 1000;

const state = {
  manifest: null,
  phenotypes: [],
  phenotypeById: new Map(),
  genes: [],
  geneById: new Map(),
  symbolIndex: {},
  topHits: [],
  qc: [],
  siteContent: null,
  currentRows: [],
  currentPairs: [],
  currentMetaPairs: [],
  currentCohortPairs: [],
  currentVisibleCohortPairs: [],
  cohortFocus: null,
  selectedPairKey: "",
  pendingFocus: null,
  lastLoadedRows: [],
  chunkCache: new Map(),
  sort: { key: "best_p", direction: "asc" },
};

const els = {};

function $(id) {
  return document.getElementById(id);
}

function formatInt(value) {
  const number = Number(value);
  return Number.isFinite(number) ? number.toLocaleString() : "NA";
}

function formatP(value) {
  const number = Number(value);
  if (!Number.isFinite(number)) return "NA";
  if (number === 0) return "0";
  if (number < 0.001) return number.toExponential(2);
  return number.toPrecision(3);
}

function formatNumber(value) {
  const number = Number(value);
  if (!Number.isFinite(number)) return "NA";
  if (Math.abs(number) >= 1000) return number.toLocaleString(undefined, { maximumFractionDigits: 1 });
  if (Math.abs(number) < 0.001 && number !== 0) return number.toExponential(2);
  return number.toPrecision(4);
}

function normalizeText(value) {
  return String(value || "").trim();
}

function isMissing(value) {
  const text = normalizeText(value);
  return !text || text === "NA" || text === "N/A" || text.toLowerCase() === "nan";
}

function sourceLabel(source) {
  if (source === "meta") return "Meta";
  if (source === "pre_meta") return "Cohort";
  return source || "NA";
}

function get(row, keys) {
  for (const key of keys) {
    if (row[key] !== undefined && row[key] !== null && row[key] !== "") return row[key];
  }
  return "";
}

function rowSource(row) {
  return get(row, ["SOURCE", "source"]);
}

function rowGeneId(row) {
  return get(row, ["ID", "gene_id"]);
}

function rowSymbol(row) {
  return get(row, ["SYMBOL", "symbol"]);
}

function rowTrait(row) {
  return get(row, ["TRAIT", "trait"]);
}

function rowP(row) {
  const value = get(row, ["P", "p"]);
  if (isMissing(value)) return Number.NaN;
  const number = Number(value);
  return Number.isFinite(number) ? number : Number.NaN;
}

function modelCode(row) {
  const code = normalizeText(get(row, ["ENCODING", "encoding"])).toUpperCase();
  if (code === "A") return "add";
  if (code === "R") return "rec";
  return code ? code.toLowerCase() : "unknown";
}

function modelLabel(code) {
  if (code === "A" || code === "add") return "Additive";
  if (code === "R" || code === "rec") return "Recessive";
  if (code === "combined") return "Combined";
  return code || "NA";
}

function modelShort(code) {
  if (code === "A" || code === "add") return "Add";
  if (code === "R" || code === "rec") return "Rec";
  if (code === "combined") return "Combined";
  return code || "NA";
}

function numericFromRow(row, keys) {
  if (!row) return Number.NaN;
  const value = get(row, keys);
  if (isMissing(value)) return Number.NaN;
  const number = Number(value);
  return Number.isFinite(number) ? number : Number.NaN;
}

function effectValue(row) {
  return numericFromRow(row, ["BETA", "beta", "Z", "z"]);
}

function rowN(row) {
  return numericFromRow(row, ["N", "n", "N_EFF", "n_eff", "number_of_pvals"]);
}

function traitLabel(id) {
  const record = state.phenotypeById.get(id);
  return record ? record.label : id;
}

function escapeHtml(value) {
  return String(value ?? "")
    .replaceAll("&", "&amp;")
    .replaceAll("<", "&lt;")
    .replaceAll(">", "&gt;")
    .replaceAll('"', "&quot;");
}

function padBucket(bucket) {
  return String(bucket).padStart(3, "0");
}

async function fetchJson(path) {
  const response = await fetch(path);
  if (!response.ok) throw new Error(`Could not load ${path}`);
  return response.json();
}

async function fetchTextMaybeGzip(path) {
  const response = await fetch(path);
  if (!response.ok) throw new Error(`Could not load ${path}`);
  const buffer = await response.arrayBuffer();
  const bytes = new Uint8Array(buffer);
  const isGzip = bytes.length >= 2 && bytes[0] === 0x1f && bytes[1] === 0x8b;
  if (!isGzip) return new TextDecoder().decode(buffer);
  if (!("DecompressionStream" in window)) {
    throw new Error("This browser cannot decompress gzip chunks.");
  }
  const stream = new Blob([buffer]).stream().pipeThrough(new DecompressionStream("gzip"));
  return new Response(stream).text();
}

function parseTsv(text) {
  const lines = text.split(/\r?\n/).filter(Boolean);
  if (!lines.length) return [];
  const header = lines[0].split("\t");
  const rows = [];
  for (let i = 1; i < lines.length; i += 1) {
    const cells = lines[i].split("\t");
    const row = {};
    for (let j = 0; j < header.length; j += 1) {
      row[header[j]] = cells[j] ?? "";
    }
    rows.push(row);
  }
  return rows;
}

function chunkAvailable(source, bucket) {
  const chunks = state.manifest?.chunks?.[source] || [];
  return chunks.some((item) => Number(item.bucket) === Number(bucket));
}

async function loadBucket(source, bucket) {
  const key = `${source}:${bucket}`;
  if (state.chunkCache.has(key)) return state.chunkCache.get(key);
  const path = `${DATA_ROOT}/chunks/${source}/gene-buckets/${padBucket(bucket)}.tsv.gz`;
  const rows = parseTsv(await fetchTextMaybeGzip(path));
  for (const row of rows) {
    row.SOURCE = source;
    const gene = state.geneById.get(row.ID);
    if (gene && gene.symbol) row.SYMBOL = gene.symbol;
  }
  state.chunkCache.set(key, rows);
  return rows;
}

function populateSelect(select, values, labelMap = {}) {
  select.innerHTML = "";
  const all = document.createElement("option");
  all.value = "";
  all.textContent = "All";
  select.appendChild(all);
  for (const value of values) {
    const option = document.createElement("option");
    option.value = value;
    option.textContent = labelMap[value] || value;
    select.appendChild(option);
  }
}

function populateCheckboxGroup(container, values, labelMap = {}) {
  container.innerHTML = "";
  for (const value of values) {
    const label = document.createElement("label");
    label.className = "checkbox-chip";
    const input = document.createElement("input");
    input.type = "checkbox";
    input.value = value;
    input.name = container.id;
    const span = document.createElement("span");
    span.textContent = labelMap[value] || value;
    label.append(input, span);
    container.appendChild(label);
  }
}

function selectedCheckboxValues(container) {
  return Array.from(container.querySelectorAll('input[type="checkbox"]:checked')).map((input) => input.value);
}

function setCheckboxValues(container, values) {
  const wanted = new Set(values.filter(Boolean));
  container.querySelectorAll('input[type="checkbox"]').forEach((input) => {
    input.checked = wanted.has(input.value);
  });
}

function categoryValues(field) {
  const values = new Set();
  const sources = state.manifest?.sources || {};
  for (const source of Object.values(sources)) {
    const categories = source.categories || {};
    for (const value of Object.keys(categories[field] || {})) {
      if (value !== "META") values.add(value);
    }
  }
  return Array.from(values).sort();
}

function setupFilters() {
  const traitIds = state.phenotypes.map((item) => item.id).sort();
  const traitLabels = Object.fromEntries(
    state.phenotypes.map((item) => [item.id, `${item.id} - ${item.label}`]),
  );
  populateSelect(els.traitFilter, traitIds, traitLabels);
  populateSelect(els.cohortFilter, categoryValues("BIOBANK"));
  populateSelect(els.ancestryFilter, categoryValues("ANCESTRY"));
  populateCheckboxGroup(els.annotationFilter, categoryValues("ANNOTATION"));
  const modelLabels = { A: "Additive", R: "Recessive" };
  populateSelect(els.encodingFilter, categoryValues("ENCODING").filter((value) => modelLabels[value]), modelLabels);
}

function renderSummary() {
  const sources = state.manifest?.sources || {};
  const rows = Object.values(sources).reduce((sum, source) => sum + Number(source.rows || 0), 0);
  const metrics = [
    ["Rows indexed", formatInt(rows)],
    ["Phenotypes", formatInt(state.manifest?.counts?.phenotypes)],
    ["Genes", formatInt(state.manifest?.counts?.genes)],
  ];
  els.summaryGrid.innerHTML = metrics
    .map(([label, value]) => `<div class="metric"><b>${label}</b><span>${value}</span></div>`)
    .join("");
}

function renderQc() {
  if (!state.qc.length) {
    els.qcSummary.innerHTML = `<div class="empty">Run the data build to populate QC.</div>`;
    return;
  }
  els.qcSummary.innerHTML = state.qc
    .map(
      (qc) => `
        <div class="mini-row"><b>${escapeHtml(sourceLabel(qc.source))}</b><span>${formatInt(qc.rows)} rows</span></div>
        <div class="mini-row"><b>Header</b><span>${qc.header_matches_expected ? "matches" : "check"}</span></div>
        <div class="mini-row"><b>No symbol</b><span>${formatInt(qc.missing_symbol)}</span></div>
        <div class="mini-row"><b>P range</b><span>${formatP(qc.pvalue?.min)} to ${formatP(qc.pvalue?.max)}</span></div>
      `,
    )
    .join("");
}

function mailtoFor(email) {
  return `mailto:${email}?subject=BRaVa%20recessive%20browser%20feedback`;
}

function applySiteContent() {
  const content = state.siteContent || {};
  if (content.landingEyebrow && els.landingEyebrow) els.landingEyebrow.textContent = content.landingEyebrow;
  if (content.landingTitle && els.landingTitle) els.landingTitle.textContent = content.landingTitle;
  if (content.landingCopy && els.landingCopy) els.landingCopy.textContent = content.landingCopy;

  if (els.paperPrefix) els.paperPrefix.textContent = content.paperPrefix || "Read our publication for more details:";
  if (els.paperCitation && content.paperUrl) {
    const prefix = escapeHtml(content.paperPrefix || "Read our publication for more details:");
    els.paperCitation.innerHTML = `<span id="paper-prefix">${prefix}</span> <a href="${escapeHtml(content.paperUrl)}" target="_blank" rel="noreferrer">${escapeHtml(content.paperLabel || content.paperUrl)}</a>`;
    els.paperPrefix = $("paper-prefix");
  }

  const email = content.contactEmail || "kalantzis@ebi.ac.uk";
  if (els.feedbackButton) els.feedbackButton.href = mailtoFor(email);

  if (content.copyright) {
    if (els.siteFootnote) els.siteFootnote.textContent = content.copyright;
    if (els.resultsFootnote) els.resultsFootnote.textContent = content.copyright;
  }
}

function selectedSources() {
  return ["meta", "pre_meta"];
}

function activeFilters() {
  return {
    trait: els.traitFilter.value,
    cohort: els.cohortFilter.value,
    ancestry: els.ancestryFilter.value,
    annotation: selectedCheckboxValues(els.annotationFilter),
    encoding: els.encodingFilter.value,
    p: els.pFilter.value ? Number(els.pFilter.value) : null,
  };
}

function rowPassesBaseFilters(row, filters) {
  const source = rowSource(row);
  const annotation = get(row, ["ANNOTATION", "annotation"]);
  if (source !== "meta" && annotation === "cauchy") return false;
  if (filters.trait && rowTrait(row) !== filters.trait) return false;
  if (source !== "meta" && filters.cohort && get(row, ["BIOBANK", "biobank"]) !== filters.cohort) return false;
  if (source !== "meta" && filters.ancestry && get(row, ["ANCESTRY", "ancestry"]) !== filters.ancestry) return false;
  if (filters.annotation.length && !filters.annotation.includes(get(row, ["ANNOTATION", "annotation"]))) return false;
  return true;
}

function pairKey(row) {
  return [
    rowSource(row),
    rowGeneId(row),
    rowTrait(row),
    get(row, ["ANNOTATION", "annotation"]),
    get(row, ["SEX", "sex"]),
    get(row, ["BIOBANK", "biobank"]),
    get(row, ["ANCESTRY", "ancestry"]),
  ].join("\u001f");
}

function seedPair(row) {
  return {
    key: pairKey(row),
    source: rowSource(row),
    geneId: rowGeneId(row),
    symbol: rowSymbol(row),
    trait: rowTrait(row),
    annotation: get(row, ["ANNOTATION", "annotation"]),
    sex: get(row, ["SEX", "sex"]),
    biobank: get(row, ["BIOBANK", "biobank"]),
    ancestry: get(row, ["ANCESTRY", "ancestry"]),
    add: null,
    rec: null,
    combined: null,
    rawRows: [],
  };
}

function setPairModelRow(pair, row) {
  const code = modelCode(row);
  const slot = code === "add" || code === "rec" || code === "combined" ? code : "combined";
  const existing = pair[slot];
  if (!existing || rowP(row) < rowP(existing)) pair[slot] = row;
  if (!pair.symbol && rowSymbol(row)) pair.symbol = rowSymbol(row);
  pair.rawRows.push(row);
}

function pairRows(rows) {
  const pairs = new Map();
  for (const row of rows) {
    const key = pairKey(row);
    if (!pairs.has(key)) pairs.set(key, seedPair(row));
    setPairModelRow(pairs.get(key), row);
  }
  return Array.from(pairs.values());
}

function pairModelRows(pair) {
  return [pair.add, pair.rec, pair.combined].filter(Boolean);
}

function pairBestP(pair) {
  const values = pairModelRows(pair).map(rowP).filter(Number.isFinite);
  return values.length ? Math.min(...values) : Number.NaN;
}

function pairPrimaryRow(pair) {
  const rows = pairModelRows(pair);
  if (!rows.length) return null;
  return rows.slice().sort((a, b) => {
    const ap = rowP(a);
    const bp = rowP(b);
    if (!Number.isFinite(ap) && !Number.isFinite(bp)) return 0;
    if (!Number.isFinite(ap)) return 1;
    if (!Number.isFinite(bp)) return -1;
    return ap - bp;
  })[0];
}

function pairN(pair) {
  const values = pairModelRows(pair).map(rowN).filter(Number.isFinite);
  return values.length ? Math.max(...values) : Number.NaN;
}

function pairCases(pair) {
  const values = pairModelRows(pair).map((row) => numericFromRow(row, ["N_CASE", "n_case"])).filter(Number.isFinite);
  return values.length ? Math.max(...values) : Number.NaN;
}

function pairControls(pair) {
  const values = pairModelRows(pair).map((row) => numericFromRow(row, ["N_CTRL", "n_ctrl"])).filter(Number.isFinite);
  return values.length ? Math.max(...values) : Number.NaN;
}

function pairDelta(pair) {
  const add = effectValue(pair.add);
  const rec = effectValue(pair.rec);
  return Number.isFinite(add) && Number.isFinite(rec) ? rec - add : Number.NaN;
}

function pairIsCauchy(pair) {
  return normalizeText(pair.annotation).toLowerCase() === "cauchy";
}

function pairNCohorts(pair) {
  return pairIsCauchy(pair) ? Number.NaN : pairN(pair);
}

function pairEnrichment(pair) {
  const recP = rowP(pair.rec || {});
  const addP = rowP(pair.add || {});
  if (!Number.isFinite(recP) || !Number.isFinite(addP) || recP <= 0 || addP <= 0) return Number.NaN;
  const recLog = -Math.log10(recP);
  const addLog = -Math.log10(addP);
  return addLog > 0 ? recLog / addLog : Number.NaN;
}

function pairPassesModelFilter(pair, filters) {
  const threshold = filters.p;
  if (filters.encoding) {
    const wanted = filters.encoding === "A" ? pair.add : pair.rec;
    if (!wanted) return false;
    if (threshold !== null) {
      const p = rowP(wanted);
      return Number.isFinite(p) && p <= threshold;
    }
    return true;
  }
  if (threshold !== null) {
    const p = pairBestP(pair);
    return Number.isFinite(p) && p <= threshold;
  }
  return true;
}

function pairSortValue(pair, key) {
  switch (key) {
    case "gene":
      return pair.symbol || pair.geneId;
    case "gene_id":
      return pair.geneId;
    case "trait":
      return traitLabel(pair.trait);
    case "source":
      return sourceLabel(pair.source);
    case "cohort":
      return pair.biobank || "META";
    case "ancestry":
      return pair.ancestry || "META";
    case "annotation":
      return pair.annotation;
    case "sex":
      return pair.sex;
    case "add_effect":
      return effectValue(pair.add);
    case "add_se":
      return numericFromRow(pair.add, ["SE", "se"]);
    case "add_p":
      return rowP(pair.add || {});
    case "rec_effect":
      return effectValue(pair.rec);
    case "rec_se":
      return numericFromRow(pair.rec, ["SE", "se"]);
    case "rec_p":
      return rowP(pair.rec || {});
    case "delta":
      return pairDelta(pair);
    case "best_p":
      return pairBestP(pair);
    case "enrichment":
      return pairEnrichment(pair);
    case "n_cohorts":
      return pairNCohorts(pair);
    case "n":
      return pairN(pair);
    case "cases":
      return pairCases(pair);
    case "controls":
      return pairControls(pair);
    default:
      return pairBestP(pair);
  }
}

function compareSortValues(left, right, direction) {
  const leftMissing = left === null || left === undefined || left === "" || (typeof left === "number" && !Number.isFinite(left));
  const rightMissing = right === null || right === undefined || right === "" || (typeof right === "number" && !Number.isFinite(right));
  if (leftMissing && rightMissing) return 0;
  if (leftMissing) return 1;
  if (rightMissing) return -1;

  let comparison;
  if (typeof left === "number" && typeof right === "number") comparison = left - right;
  else comparison = String(left).localeCompare(String(right), undefined, { numeric: true, sensitivity: "base" });
  return direction === "desc" ? -comparison : comparison;
}

function sortPairs(pairs) {
  const { key, direction } = state.sort;
  return pairs.slice().sort((a, b) => {
    const primary = compareSortValues(pairSortValue(a, key), pairSortValue(b, key), direction);
    if (primary !== 0) return primary;
    return compareSortValues(pairBestP(a), pairBestP(b), "asc");
  });
}

function applyFilters(rows) {
  const filters = activeFilters();
  const baseRows = rows.filter((row) => rowPassesBaseFilters(row, filters));
  return sortPairs(pairRows(baseRows).filter((pair) => pairPassesModelFilter(pair, filters)));
}

function updateUrl(query = "", options = {}) {
  const params = new URLSearchParams();
  const includeFocus = options.focus !== false;
  if (query) params.set("q", query);
  for (const [key, value] of Object.entries(activeFilters())) {
    if (Array.isArray(value)) {
      if (value.length) params.set(key, value.join(","));
    } else if (value !== "" && value !== null) {
      params.set(key, value);
    }
  }
  const focus = state.cohortFocus || state.pendingFocus;
  if (includeFocus && focus?.geneId && focus?.trait) {
    params.set("focus_gene", focus.geneId);
    params.set("focus_trait", focus.trait);
  }
  const next = `${window.location.pathname}${params.toString() ? `?${params.toString()}` : ""}`;
  window.history.replaceState({}, "", next);
}

function showMessage(text) {
  if (!text) {
    els.message.hidden = true;
    els.message.textContent = "";
    return;
  }
  els.message.hidden = false;
  els.message.textContent = text;
}


function showLandingPage() {
  document.body.classList.add("home-view");
  document.body.classList.remove("results-view");
  if (els.landingPage) els.landingPage.hidden = false;
  if (els.resultsPage) els.resultsPage.hidden = true;
}

function showResultsPage() {
  document.body.classList.remove("home-view");
  document.body.classList.add("results-view");
  if (els.landingPage) els.landingPage.hidden = true;
  if (els.resultsPage) els.resultsPage.hidden = false;
}

function renderRows(rows, title, subtitle) {
  const pairs = applyFilters(rows);
  const metaPairs = pairs.filter((pair) => pair.source === "meta");
  const cohortPairs = pairs.filter((pair) => pair.source !== "meta");
  state.currentRows = rows;
  state.currentPairs = pairs;
  state.currentMetaPairs = metaPairs;
  state.currentCohortPairs = cohortPairs;
  state.cohortFocus = null;
  state.selectedPairKey = "";
  els.resultTitle.textContent = title;
  const shownMeta = Math.min(MAX_TABLE_ROWS, metaPairs.length);
  const shownCohort = Math.min(MAX_TABLE_ROWS, cohortPairs.length);
  els.resultSubtitle.textContent = `${subtitle} ${metaPairs.length.toLocaleString()} meta-analysis pairs and ${cohortPairs.length.toLocaleString()} cohort-level pairs match. Showing ${shownMeta.toLocaleString()} meta and ${shownCohort.toLocaleString()} cohort rows.`;
  renderTable("meta", metaPairs.slice(0, MAX_TABLE_ROWS));
  renderCohortPanel();
  if (!applyPendingFocus()) renderFocusedPlot();
}

function pairFocusFromPair(pair) {
  if (!pair) return null;
  return {
    geneId: pair.geneId || "",
    trait: pair.trait || "",
    symbol: pair.symbol || "",
  };
}

function focusMatchesPair(focus, pair) {
  return Boolean(focus && pair && pair.geneId === focus.geneId && pair.trait === focus.trait);
}

function focusLabel(focus) {
  if (!focus) return "all matching gene-phenotype pairs";
  const gene = focus.symbol || focus.geneId || "gene";
  return `${focus.trait || "phenotype"}-${gene}`;
}

function cohortPairsForCurrentFocus() {
  if (!state.cohortFocus) return state.currentCohortPairs;
  const filters = activeFilters();
  const rows = state.currentRows.filter(
    (row) =>
      rowSource(row) !== "meta" &&
      rowGeneId(row) === state.cohortFocus.geneId &&
      rowTrait(row) === state.cohortFocus.trait &&
      rowPassesBaseFilters(row, filters),
  );
  return sortPairs(pairRows(rows));
}


function rowMergeKey(row) {
  return [
    rowSource(row),
    rowGeneId(row),
    rowTrait(row),
    get(row, ["BIOBANK", "biobank"]),
    get(row, ["ANCESTRY", "ancestry"]),
    get(row, ["ANNOTATION", "annotation"]),
    get(row, ["SEX", "sex"]),
    get(row, ["ENCODING", "encoding"]),
  ].join("\u001f");
}

async function loadCohortRowsForPair(pair) {
  if (!pair?.geneId) return false;
  const gene = state.geneById.get(pair.geneId);
  if (!gene || gene.bucket === undefined || gene.bucket === null || !chunkAvailable("pre_meta", gene.bucket)) return false;
  const alreadyLoaded = state.currentRows.some(
    (row) => rowSource(row) === "pre_meta" && rowGeneId(row) === pair.geneId && rowTrait(row) === pair.trait,
  );
  if (alreadyLoaded) return false;

  const previousStatus = els.dataStatus ? els.dataStatus.textContent : "";
  if (els.dataStatus) els.dataStatus.textContent = `Loading cohort rows for ${pair.symbol || pair.geneId}`;
  const loadedRows = (await loadBucket("pre_meta", gene.bucket)).filter((row) => rowGeneId(row) === pair.geneId);
  if (!loadedRows.length) {
    if (els.dataStatus) els.dataStatus.textContent = previousStatus;
    return false;
  }

  const seen = new Set(state.currentRows.map(rowMergeKey));
  const addedRows = loadedRows.filter((row) => {
    const key = rowMergeKey(row);
    if (seen.has(key)) return false;
    seen.add(key);
    return true;
  });
  if (!addedRows.length) {
    if (els.dataStatus) els.dataStatus.textContent = previousStatus;
    return false;
  }

  state.currentRows = [...state.currentRows, ...addedRows];
  state.lastLoadedRows = state.currentRows;
  const pairs = applyFilters(state.currentRows);
  state.currentPairs = pairs;
  state.currentMetaPairs = pairs.filter((candidate) => candidate.source === "meta");
  state.currentCohortPairs = pairs.filter((candidate) => candidate.source !== "meta");
  const shownMeta = Math.min(MAX_TABLE_ROWS, state.currentMetaPairs.length);
  const shownCohort = Math.min(MAX_TABLE_ROWS, state.currentCohortPairs.length);
  els.resultSubtitle.textContent = `Loaded cohort rows for ${pair.symbol || pair.geneId}. ${state.currentMetaPairs.length.toLocaleString()} meta-analysis pairs and ${state.currentCohortPairs.length.toLocaleString()} cohort-level pairs match. Showing ${shownMeta.toLocaleString()} meta and ${shownCohort.toLocaleString()} cohort rows.`;
  if (els.dataStatus) els.dataStatus.textContent = `${addedRows.length.toLocaleString()} cohort rows added`;
  return true;
}

async function handlePairSelection(pair, kind) {
  if (!pair) return;
  renderDetail(pair);
  if (kind === "meta") await loadCohortRowsForPair(pair);
  setCohortFocus(pair, pair);
}

function renderCohortPanel() {
  const visiblePairs = cohortPairsForCurrentFocus();
  state.currentVisibleCohortPairs = visiblePairs;
  if (els.cohortFocusLabel) {
    els.cohortFocusLabel.textContent = state.cohortFocus
      ? `showing ${focusLabel(state.cohortFocus)}`
      : "showing all matching pairs";
  }
  if (els.clearCohortFocusButton) {
    els.clearCohortFocusButton.hidden = !state.cohortFocus;
  }
  renderTable("cohort", visiblePairs.slice(0, MAX_TABLE_ROWS));
}

function setCohortFocus(pair, selectedPair = pair, options = {}) {
  state.cohortFocus = pairFocusFromPair(pair);
  state.selectedPairKey = selectedPair?.key || pair?.key || "";
  renderTable("meta", state.currentMetaPairs.slice(0, MAX_TABLE_ROWS));
  renderCohortPanel();
  renderFocusedPlot(selectedPair);
  if (options.updateUrl !== false) updateUrl(els.searchInput.value.trim());
}

function clearCohortFocus() {
  state.cohortFocus = null;
  state.selectedPairKey = "";
  renderTable("meta", state.currentMetaPairs.slice(0, MAX_TABLE_ROWS));
  renderCohortPanel();
  renderFocusedPlot();
  updateUrl(els.searchInput.value.trim(), { focus: false });
}

function applyPendingFocus() {
  if (!state.pendingFocus) return false;
  const focus = state.pendingFocus;
  const pair =
    state.currentMetaPairs.find((candidate) => candidate.geneId === focus.geneId && candidate.trait === focus.trait) ||
    state.currentCohortPairs.find((candidate) => candidate.geneId === focus.geneId && candidate.trait === focus.trait);
  state.pendingFocus = null;
  if (!pair) return false;
  setCohortFocus(pair, pair, { updateUrl: false });
  return true;
}

function renderFocusedPlot(selected = null) {
  if (els.plotTitle) {
    els.plotTitle.textContent = state.cohortFocus
      ? `Recessive effect estimates across cohorts for ${focusLabel(state.cohortFocus)}`
      : "Recessive effect estimates across cohorts";
  }
  if (!state.cohortFocus) {
    if (els.plotPanel) {
      els.plotPanel.classList.add("empty");
      els.plotPanel.textContent = "Select a meta-analysis row or cohort-level row to plot one gene-phenotype pair.";
    }
    return;
  }
  renderPlot(cohortPairsForCurrentFocus(), selected);
}

function tableColumns(kind) {
  if (kind === "meta") {
    return [
      { key: "gene", label: "Gene symbol" },
      { key: "gene_id", label: "Gene ID" },
      { key: "trait", label: "Phenotype label" },
      { key: "annotation", label: "Annotation" },
      { key: "sex", label: "Sex" },
      { key: "rec_p", label: "Rec P" },
      { key: "add_p", label: "Add P" },
      { key: "enrichment", label: "Enrichment" },
      { key: "n_cohorts", label: "N cohorts" },
    ];
  }
  return [
    { key: "gene", label: "Gene symbol" },
    { key: "gene_id", label: "Gene ID" },
    { key: "trait", label: "Phenotype label" },
    { key: "cohort", label: "Cohort" },
    { key: "ancestry", label: "Ancestry" },
    { key: "annotation", label: "Annotation" },
    { key: "sex", label: "Sex" },
    { key: "rec_p", label: "Rec P" },
    { key: "rec_effect", label: "Rec BETA/Z" },
    { key: "rec_se", label: "Rec SE" },
    { key: "add_p", label: "Add P" },
    { key: "add_effect", label: "Add BETA/Z" },
    { key: "add_se", label: "Add SE" },
    { key: "delta", label: "Rec-Add" },
    { key: "n", label: "N" },
    { key: "cases", label: "Cases" },
    { key: "controls", label: "Controls" },
  ];
}

function sortIndicator(key) {
  if (state.sort.key !== key) return "";
  return state.sort.direction === "asc" ? " ↑" : " ↓";
}

function renderModelMetric(row, metric) {
  if (!row) return '<span class="empty">NA</span>';
  if (metric === "effect") return escapeHtml(formatNumber(effectValue(row)));
  if (metric === "se") return escapeHtml(formatNumber(numericFromRow(row, ["SE", "se"])));
  if (metric === "p") return escapeHtml(formatP(rowP(row)));
  return "NA";
}

function tableCell(pair, column) {
  const geneLabel = pair.symbol || "NA";
  switch (column.key) {
    case "gene":
      return `<b>${escapeHtml(geneLabel)}</b>`;
    case "gene_id":
      return escapeHtml(pair.geneId || "NA");
    case "trait":
      return escapeHtml(traitLabel(pair.trait));
    case "cohort":
      return escapeHtml(pair.biobank || "NA");
    case "ancestry":
      return escapeHtml(pair.ancestry || "NA");
    case "annotation":
      return escapeHtml(pair.annotation || "NA");
    case "sex":
      return escapeHtml(pair.sex || "NA");
    case "rec_p":
      return `<span class="model-cell-inline rec">${renderModelMetric(pair.rec, "p")}</span>`;
    case "rec_effect":
      return renderModelMetric(pair.rec, "effect");
    case "rec_se":
      return renderModelMetric(pair.rec, "se");
    case "add_p":
      return `<span class="model-cell-inline add">${renderModelMetric(pair.add, "p")}</span>`;
    case "add_effect":
      return renderModelMetric(pair.add, "effect");
    case "add_se":
      return renderModelMetric(pair.add, "se");
    case "delta":
      return escapeHtml(formatNumber(pairDelta(pair)));
    case "best_p":
      return escapeHtml(formatP(pairBestP(pair)));
    case "enrichment":
      return escapeHtml(formatNumber(pairEnrichment(pair)));
    case "n_cohorts":
      return escapeHtml(formatInt(pairNCohorts(pair)));
    case "n":
      return escapeHtml(formatInt(pairN(pair)));
    case "cases":
      return escapeHtml(formatInt(pairCases(pair)));
    case "controls":
      return escapeHtml(formatInt(pairControls(pair)));
    default:
      return "NA";
  }
}

function renderTable(kind, pairs) {
  const table = kind === "meta" ? els.metaTable : els.cohortTable;
  const head = table ? table.querySelector("thead") : null;
  const body = table ? table.querySelector("tbody") : null;
  if (!head || !body) return;
  const columns = tableColumns(kind);
  head.innerHTML = `<tr>${columns
    .map(
      (column) =>
        `<th><button type="button" class="sort-header" data-sort-kind="${kind}" data-sort-key="${column.key}" aria-label="Sort by ${escapeHtml(column.label)}">${escapeHtml(column.label)}${sortIndicator(column.key)}</button></th>`,
    )
    .join("")}</tr>`;
  head.querySelectorAll("button[data-sort-key]").forEach((button) => {
    button.addEventListener("click", () => {
      const key = button.dataset.sortKey;
      if (state.sort.key === key) state.sort.direction = state.sort.direction === "asc" ? "desc" : "asc";
      else state.sort = { key, direction: key.endsWith("_p") || key === "best_p" ? "asc" : "asc" };
      if (button.dataset.sortKind === "meta") {
        state.currentMetaPairs = sortPairs(state.currentMetaPairs);
        renderTable("meta", state.currentMetaPairs.slice(0, MAX_TABLE_ROWS));
      } else {
        state.currentCohortPairs = sortPairs(state.currentCohortPairs);
        renderCohortPanel();
        renderFocusedPlot();
      }
    });
  });

  if (!pairs.length) {
    body.innerHTML = `<tr><td colspan="${columns.length}" class="empty">No ${kind === "meta" ? "meta-analysis" : "cohort-level"} association pairs match the current filters.</td></tr>`;
    return;
  }

  body.innerHTML = pairs
    .map((pair, index) => {
      const isSelected = state.selectedPairKey && pair.key === state.selectedPairKey;
      const isFocusedCohort = kind === "cohort" && state.cohortFocus && focusMatchesPair(state.cohortFocus, pair);
      const rowClass = isSelected || isFocusedCohort ? ' class="selected-row"' : "";
      return `
        <tr${rowClass} data-row-kind="${kind}" data-row-index="${index}">
          ${columns.map((column) => `<td>${tableCell(pair, column)}</td>`).join("")}
        </tr>
      `;
    })
    .join("");
  body.querySelectorAll("tr[data-row-index]").forEach((rowEl) => {
    rowEl.addEventListener("click", () => {
      const kind = rowEl.dataset.rowKind;
      const pairList = kind === "meta" ? state.currentMetaPairs : state.currentVisibleCohortPairs;
      const pair = pairList[Number(rowEl.dataset.rowIndex)];
      handlePairSelection(pair, kind);
    });
  });
}

function renderDetail(pair) {
  const phenotype = state.phenotypeById.get(pair.trait);
  const primary = pairPrimaryRow(pair) || {};
  const entries = [
    ["Gene", pair.symbol ? `${pair.symbol} (${pair.geneId})` : pair.geneId],
    ["Phenotype", phenotype ? `${pair.trait} - ${phenotype.label}` : pair.trait],
    ["Type", phenotype?.type || "NA"],
    ["Source", sourceLabel(pair.source)],
    ["Cohort", pair.biobank || "META"],
    ["Ancestry", pair.ancestry || "META"],
    ["Annotation", pair.annotation || "NA"],
    ["Sex", pair.sex || "NA"],
    ["Recessive P", formatP(rowP(pair.rec || {}))],
    ["Recessive BETA/Z", formatNumber(effectValue(pair.rec))],
    ["Recessive SE", formatNumber(numericFromRow(pair.rec, ["SE", "se"]))],
    ["Additive P", formatP(rowP(pair.add || {}))],
    ["Additive BETA/Z", formatNumber(effectValue(pair.add))],
    ["Additive SE", formatNumber(numericFromRow(pair.add, ["SE", "se"]))],
    ["Rec-Add", formatNumber(pairDelta(pair))],
    ["Best P", formatP(pairBestP(pair))],
    ["Combined P", pair.combined ? formatP(rowP(pair.combined)) : "NA"],
    ["N", formatInt(pairN(pair))],
    ["Cases", formatInt(pairCases(pair))],
    ["Controls", formatInt(pairControls(pair))],
    ["AF", formatNumber(numericFromRow(primary, ["AF"]))],
    ["AC", formatInt(numericFromRow(primary, ["AC"]))],
  ];
  els.detailPanel.classList.remove("empty");
  els.detailPanel.innerHTML = entries
    .map(([label, value]) => `<div class="mini-row"><b>${escapeHtml(label)}</b><span>${escapeHtml(value)}</span></div>`)
    .join("");
}

function ancestryOrder(value) {
  const order = ["EUR", "SAS", "AFR", "AMR", "EAS", "META"];
  const index = order.indexOf(value);
  return index === -1 ? order.length : index;
}

const PLOT_BIOBANKS = [
  { keys: ["UKB"], label: "UKB", color: "#1f77b4" },
  { keys: ["AOU"], label: "AOU", color: "#ff7f0e" },
  { keys: ["GEL", "100KGP", "100K"], label: "100kGP", color: "#2ca02c" },
  { keys: ["GNH", "G&H", "G AND H"], label: "G&H", color: "#e41a1c" },
  { keys: ["BBJ"], label: "BBJ", color: "#9467bd" },
  { keys: ["BIOME", "BIO ME"], label: "BioMe", color: "#8c564b" },
];

const PLOT_ANNOTATIONS = [
  { key: "nonsynonymous", label: "nonsynonymous", shape: "triangle" },
  { key: "pLoF", label: "pLoF", shape: "square" },
  { key: "pLoF_damaging_missense", label: "pLoF+damaging missense/protein altering variants", shape: "circle" },
];

function plotBiobankValue(pair) {
  return normalizeText(pair.biobank || (pair.source === "meta" ? "META" : "NA"));
}

function plotBiobankSpec(value) {
  const raw = normalizeText(value);
  const upper = raw.toUpperCase();
  const match = PLOT_BIOBANKS.find((entry) => entry.keys.some((key) => upper === key || upper.includes(key)));
  if (match) return match;
  if (upper === "META") return { keys: ["META"], label: "Meta", color: "#111827" };
  return { keys: [upper || "NA"], label: raw || "NA", color: "#64748b" };
}

function plotBiobankOrder(value) {
  const spec = plotBiobankSpec(value);
  const index = PLOT_BIOBANKS.findIndex((entry) => entry.label === spec.label);
  return index === -1 ? PLOT_BIOBANKS.length : index;
}

function plotAnnotationSpec(value) {
  const raw = normalizeText(value);
  const lower = raw.toLowerCase();
  if (lower === "nonsynonymous") return PLOT_ANNOTATIONS[0];
  if (lower === "plof") return PLOT_ANNOTATIONS[1];
  if (lower.includes("damaging") || lower.includes("protein")) return PLOT_ANNOTATIONS[2];
  return { key: lower || "other", label: raw || "Other", shape: "circle" };
}

function plotAnnotationOrder(value) {
  const spec = plotAnnotationSpec(value);
  const index = PLOT_ANNOTATIONS.findIndex((entry) => entry.key === spec.key);
  return index === -1 ? PLOT_ANNOTATIONS.length : index;
}

function plotMarkerShape(shape, x, y, color, selected, size = 6) {
  const stroke = selected ? "#111827" : color;
  const strokeWidth = selected ? 1.5 : 0;
  if (shape === "triangle") {
    return `<polygon points="${x},${y - size} ${x - size},${y + size} ${x + size},${y + size}" fill="${color}" stroke="${stroke}" stroke-width="${strokeWidth}" />`;
  }
  if (shape === "square") {
    const side = size * 1.65;
    return `<rect x="${x - side / 2}" y="${y - side / 2}" width="${side}" height="${side}" fill="${color}" stroke="${stroke}" stroke-width="${strokeWidth}" />`;
  }
  return `<circle cx="${x}" cy="${y}" r="${selected ? size * 0.9 : size * 0.78}" fill="${color}" stroke="${stroke}" stroke-width="${strokeWidth}" />`;
}

function plotPointMarker(point, x, y, selected) {
  return plotMarkerShape(point.annotationSpec.shape, x, y, point.biobankSpec.color, selected);
}

function plotLineStyle(point) {
  const p = rowP(point.row);
  return Number.isFinite(p) && p < 0.05 ? "" : 'stroke-dasharray="4 4"';
}

function plotLegendSvg() {
  const biobankRows = [
    [PLOT_BIOBANKS[0], PLOT_BIOBANKS[3]],
    [PLOT_BIOBANKS[1], PLOT_BIOBANKS[4]],
    [PLOT_BIOBANKS[2], PLOT_BIOBANKS[5]],
  ];
  const biobankItems = biobankRows
    .map((row, index) =>
      row
        .map((entry, column) => {
          const x = 102 + column * 92;
          const y = 42 + index * 22;
          return `<circle cx="${x}" cy="${y - 4}" r="6" fill="${entry.color}" /><text x="${x + 18}" y="${y}" font-size="12" fill="#344054">${escapeHtml(entry.label)}</text>`;
        })
        .join(""),
    )
    .join("");

  const annotationItems = PLOT_ANNOTATIONS.map((entry, index) => {
    const y = 42 + index * 24;
    const label = entry.key === "pLoF_damaging_missense"
      ? `<text x="318" y="${y - 2}" font-size="12" fill="#344054">pLoF+damaging missense/</text><text x="318" y="${y + 12}" font-size="12" fill="#344054">protein altering variants</text>`
      : `<text x="318" y="${y}" font-size="12" fill="#344054">${escapeHtml(entry.label)}</text>`;
    return `${plotMarkerShape(entry.shape, 294, y - 4, "#111827", false)}${label}`;
  }).join("");

  return `
    <g aria-label="Forest plot legend">
      <text x="76" y="18" font-size="14" font-weight="800" fill="#111827">Biobank</text>
      ${biobankItems}
      <text x="274" y="18" font-size="14" font-weight="800" fill="#111827">Annotation</text>
      ${annotationItems}
      <text x="600" y="18" font-size="14" font-weight="800" fill="#111827">P &lt; 0.05</text>
      <line x1="600" x2="650" y1="42" y2="42" stroke="#111827" stroke-width="2" />
      <text x="662" y="46" font-size="12" fill="#344054">Yes</text>
      <line x1="600" x2="650" y1="66" y2="66" stroke="#111827" stroke-width="2" stroke-dasharray="4 4" />
      <text x="662" y="70" font-size="12" fill="#344054">No</text>
    </g>
  `;
}

function renderPlot(pairs, selected = null) {
  const points = [];
  for (const pair of pairs) {
    const row = pair.rec;
    if (!row) continue;
    const beta = effectValue(row);
    const se = numericFromRow(row, ["SE", "se"]);
    if (!Number.isFinite(beta) || !Number.isFinite(se)) continue;
    const ancestry = pair.ancestry || (pair.source === "meta" ? "META" : "NA");
    const biobank = plotBiobankValue(pair);
    const annotation = pair.annotation || get(row, ["ANNOTATION", "annotation"]);
    points.push({
      pair,
      row,
      beta,
      se,
      lo: beta - 1.96 * se,
      hi: beta + 1.96 * se,
      ancestry,
      biobank,
      annotation,
      biobankSpec: plotBiobankSpec(biobank),
      annotationSpec: plotAnnotationSpec(annotation),
    });
    if (points.length >= 84) break;
  }

  if (!points.length) {
    els.plotPanel.classList.add("empty");
    els.plotPanel.textContent = "Rows with recessive beta and SE will be plotted here.";
    return;
  }

  const grouped = new Map();
  for (const point of points) {
    if (!grouped.has(point.ancestry)) grouped.set(point.ancestry, []);
    grouped.get(point.ancestry).push(point);
  }
  const groups = Array.from(grouped, ([ancestry, values]) => ({ ancestry, values }))
    .sort((a, b) => ancestryOrder(a.ancestry) - ancestryOrder(b.ancestry) || a.ancestry.localeCompare(b.ancestry));

  const width = 780;
  const rowHeight = 18;
  const groupPadding = 18;
  const groupGap = 6;
  const legendHeight = 96;
  const margin = { top: legendHeight + 18, right: 30, bottom: 38, left: 94 };
  const groupHeights = groups.map((group) => Math.max(50, group.values.length * rowHeight + groupPadding * 2));
  const plotHeight = groupHeights.reduce((sum, height) => sum + height, 0) + Math.max(0, groups.length - 1) * groupGap;
  const height = margin.top + margin.bottom + plotHeight;
  let minX = Math.min(-0.01, ...points.map((point) => point.lo));
  let maxX = Math.max(0.01, ...points.map((point) => point.hi));
  if (minX === maxX) {
    minX -= 0.1;
    maxX += 0.1;
  } else {
    const padding = (maxX - minX) * 0.08;
    minX -= padding;
    maxX += padding;
  }
  const scale = (value) => margin.left + ((value - minX) / (maxX - minX)) * (width - margin.left - margin.right);
  const zero = scale(0);
  let yCursor = margin.top;

  const layers = groups
    .map((group, groupIndex) => {
      const groupHeight = groupHeights[groupIndex];
      const band = groupIndex % 2 === 0 ? "#f3f4f6" : "#ffffff";
      const yStart = yCursor;
      const yMid = yStart + groupHeight / 2;
      const items = group.values
        .sort((a, b) => {
          const cohortCompare = plotBiobankOrder(a.biobank) - plotBiobankOrder(b.biobank);
          if (cohortCompare !== 0) return cohortCompare;
          const annotationCompare = plotAnnotationOrder(a.annotation) - plotAnnotationOrder(b.annotation);
          if (annotationCompare !== 0) return annotationCompare;
          return compareSortValues(rowP(a.row), rowP(b.row), "asc");
        })
        .map((point, idx) => {
          const y = yStart + groupPadding + idx * rowHeight + rowHeight / 2;
          const isSelected = selected && selected.key === point.pair.key;
          const color = point.biobankSpec.color;
          const label = `${point.biobankSpec.label} ${point.annotationSpec.label}`;
          return `
            <title>${escapeHtml(group.ancestry)} ${escapeHtml(label)} recessive beta ${formatNumber(point.beta)}, SE ${formatNumber(point.se)}, P ${formatP(rowP(point.row))}</title>
            <line x1="${scale(point.lo)}" x2="${scale(point.hi)}" y1="${y}" y2="${y}" stroke="${color}" stroke-width="${isSelected ? 2.2 : 1.6}" ${plotLineStyle(point)} opacity="0.9" />
            ${plotPointMarker(point, scale(point.beta), y, isSelected)}
          `;
        })
        .join("");
      yCursor += groupHeight + groupGap;
      return `
        <rect x="${margin.left}" y="${yStart}" width="${width - margin.left - margin.right}" height="${groupHeight}" fill="${band}" />
        <text x="10" y="${yMid + 5}" font-size="18" fill="#475467" font-weight="700">${escapeHtml(group.ancestry)}</text>
        ${items}
      `;
    })
    .join("");

  els.plotPanel.classList.remove("empty");
  els.plotPanel.innerHTML = `
    <svg viewBox="0 0 ${width} ${height}" role="img" aria-label="Recessive effect estimates across cohorts">
      ${plotLegendSvg()}
      ${layers}
      <line x1="${zero}" x2="${zero}" y1="${margin.top - 8}" y2="${height - margin.bottom + 8}" stroke="#111827" stroke-dasharray="5 5" opacity="0.65" />
      <text x="${margin.left}" y="${height - 12}" font-size="12" fill="#667085">${formatNumber(minX)}</text>
      <text x="${zero - 4}" y="${height - 12}" font-size="12" fill="#667085">0</text>
      <text x="${width - margin.right - 42}" y="${height - 12}" font-size="12" fill="#667085">${formatNumber(maxX)}</text>
      <text x="${width / 2 - 18}" y="${height - 8}" font-size="13" font-weight="700" fill="#111827">Beta</text>
    </svg>
  `;
}

function geneMatches(query) {
  const raw = query.trim();
  if (!raw) return [];
  const upper = raw.toUpperCase();
  if (state.symbolIndex[upper]) {
    return state.symbolIndex[upper].map((id) => state.geneById.get(id)).filter(Boolean);
  }
  if (state.geneById.has(raw)) return [state.geneById.get(raw)];
  if (state.geneById.has(upper)) return [state.geneById.get(upper)];
  const starts = state.genes
    .filter((gene) => gene.symbol && gene.symbol.toUpperCase().startsWith(upper))
    .slice(0, 12);
  return starts;
}

function phenotypeMatches(query) {
  const upper = query.trim().toUpperCase();
  if (!upper) return [];
  return state.phenotypes.filter(
    (item) => item.id.toUpperCase() === upper || item.label.toUpperCase().includes(upper),
  );
}

function populateSearchSuggestions() {
  if (!els.searchSuggestions) return;
  const options = [];
  for (const gene of state.genes) {
    if (gene.symbol) options.push(`<option value="${escapeHtml(gene.symbol)}" label="${escapeHtml(gene.id)}"></option>`);
    options.push(`<option value="${escapeHtml(gene.id)}" label="${escapeHtml(gene.symbol || "Ensembl gene ID")}"></option>`);
  }
  for (const phenotype of state.phenotypes) {
    options.push(`<option value="${escapeHtml(phenotype.id)}" label="${escapeHtml(phenotype.label)}"></option>`);
    if (phenotype.label) options.push(`<option value="${escapeHtml(phenotype.label)}" label="${escapeHtml(phenotype.id)}"></option>`);
  }
  els.searchSuggestions.innerHTML = options.join("");
}

async function runSearch(options = {}) {
  showResultsPage();
  const query = els.searchInput.value.trim();
  const preserveFocus = options.preserveFocus === true;
  if (!preserveFocus) state.pendingFocus = null;
  updateUrl(query, { focus: preserveFocus });
  showMessage("");
  if (!query) {
    renderRows(state.topHits, "Top associations", "Smallest available p-values from the generated index.");
    return;
  }

  const genes = geneMatches(query);
  if (genes.length) {
    els.dataStatus.textContent = "Loading gene chunk";
    const ids = new Set(genes.map((gene) => gene.id));
    const buckets = Array.from(new Set(genes.map((gene) => gene.bucket)));
    const rows = [];
    const missing = [];
    for (const source of selectedSources()) {
      for (const bucket of buckets) {
        if (!chunkAvailable(source, bucket)) {
          missing.push(`${sourceLabel(source)} bucket ${padBucket(bucket)}`);
          continue;
        }
        rows.push(...(await loadBucket(source, bucket)).filter((row) => ids.has(rowGeneId(row))));
      }
    }
    state.lastLoadedRows = rows;
    els.dataStatus.textContent = `${rows.length.toLocaleString()} rows loaded`;
    const geneLabel = genes.map((gene) => gene.symbol || gene.id).join(", ");
    if (missing.length) showMessage(`Some chunks are not present: ${missing.slice(0, 4).join(", ")}.`);
    renderRows(rows, `Gene search: ${geneLabel}`, "Rows loaded from gene-bucket chunks.");
    return;
  }

  const phenotypes = phenotypeMatches(query);
  if (phenotypes.length) {
    const ids = new Set(phenotypes.map((item) => item.id));
    els.traitFilter.value = phenotypes[0].id;
    const rows = state.topHits.filter((row) => ids.has(rowTrait(row)));
    state.lastLoadedRows = rows;
    renderRows(rows, `Phenotype search: ${phenotypes[0].id}`, "Showing indexed top associations for matching phenotype.");
    showMessage("Full phenotype-wide loading is not generated yet; this view uses the top-hit index.");
    return;
  }

  renderRows([], `Search: ${query}`, "No matching gene or phenotype was found.");
}

function renderTopHits() {
  showResultsPage();
  els.searchInput.value = "";
  state.pendingFocus = null;
  updateUrl("", { focus: false });
  showMessage("");
  state.lastLoadedRows = state.topHits;
  renderRows(state.topHits, "Top associations", "Smallest available p-values from the generated index.");
}

function resetFilters() {
  state.cohortFocus = null;
  state.selectedPairKey = "";
  state.pendingFocus = null;
  els.traitFilter.value = "";
  els.cohortFilter.value = "";
  els.ancestryFilter.value = "";
  setCheckboxValues(els.annotationFilter, []);
  els.encodingFilter.value = "";
  els.pFilter.value = "";
  renderRows(state.lastLoadedRows.length ? state.lastLoadedRows : state.topHits, els.resultTitle.textContent, "Filters reset.");
  updateUrl(els.searchInput.value.trim());
}

function pairToExportRow(pair) {
  const base = {
    gene_symbol: pair.symbol || "",
    gene_id: pair.geneId,
    phenotype_label: traitLabel(pair.trait),
    annotation: pair.annotation || "",
    sex: pair.sex || "",
    rec_p: formatP(rowP(pair.rec || {})),
    add_p: formatP(rowP(pair.add || {})),
    best_p: formatP(pairBestP(pair)),
  };
  if (pair.source === "meta") {
    return { ...base, enrichment: formatNumber(pairEnrichment(pair)), n_cohorts: formatInt(pairNCohorts(pair)) };
  }
  return {
    ...base,
    cohort: pair.biobank || "",
    ancestry: pair.ancestry || "",
    rec_effect_or_z: formatNumber(effectValue(pair.rec)),
    rec_se: formatNumber(numericFromRow(pair.rec, ["SE", "se"])),
    add_effect_or_z: formatNumber(effectValue(pair.add)),
    add_se: formatNumber(numericFromRow(pair.add, ["SE", "se"])),
    rec_minus_add: formatNumber(pairDelta(pair)),
    n: formatInt(pairN(pair)),
    cases: formatInt(pairCases(pair)),
    controls: formatInt(pairControls(pair)),
  };
}

function downloadRows() {
  const rows = [...state.currentMetaPairs, ...state.currentVisibleCohortPairs].map(pairToExportRow);
  if (!rows.length) return;
  const columns = [
    "gene_symbol",
    "gene_id",
    "phenotype_label",
    "cohort",
    "ancestry",
    "annotation",
    "sex",
    "rec_p",
    "rec_effect_or_z",
    "rec_se",
    "add_p",
    "add_effect_or_z",
    "add_se",
    "rec_minus_add",
    "best_p",
    "enrichment",
    "n_cohorts",
    "n",
    "cases",
    "controls",
  ];
  const lines = [columns.join("\t")];
  for (const row of rows) {
    lines.push(columns.map((column) => String(row[column] ?? "").replace(/\t/g, " ")).join("\t"));
  }
  const blob = new Blob([lines.join("\n") + "\n"], { type: "text/tab-separated-values" });
  const url = URL.createObjectURL(blob);
  const link = document.createElement("a");
  link.href = url;
  link.download = "brava_recessive_browser_paired_rows.tsv";
  link.click();
  URL.revokeObjectURL(url);
}

async function downloadPlotPng() {
  const svg = els.plotPanel ? els.plotPanel.querySelector("svg") : null;
  if (!svg) {
    if (els.dataStatus) els.dataStatus.textContent = "No plot to save";
    return;
  }
  const clone = svg.cloneNode(true);
  clone.setAttribute("xmlns", "http://www.w3.org/2000/svg");
  const source = new XMLSerializer().serializeToString(clone);
  const svgBlob = new Blob([source], { type: "image/svg+xml;charset=utf-8" });
  const url = URL.createObjectURL(svgBlob);
  const image = new Image();
  image.decoding = "async";
  const loaded = new Promise((resolve, reject) => {
    image.onload = resolve;
    image.onerror = reject;
  });
  image.src = url;
  try {
    await loaded;
    const box = svg.viewBox.baseVal;
    const width = Math.max(1, Math.ceil(box?.width || svg.clientWidth || 780));
    const height = Math.max(1, Math.ceil(box?.height || svg.clientHeight || 480));
    const scale = 2;
    const canvas = document.createElement("canvas");
    canvas.width = width * scale;
    canvas.height = height * scale;
    const context = canvas.getContext("2d");
    context.fillStyle = "#ffffff";
    context.fillRect(0, 0, canvas.width, canvas.height);
    context.drawImage(image, 0, 0, canvas.width, canvas.height);
    const pngUrl = canvas.toDataURL("image/png");
    const link = document.createElement("a");
    link.href = pngUrl;
    link.download = "brava_recessive_forest_plot.png";
    link.click();
  } finally {
    URL.revokeObjectURL(url);
  }
}

async function copyLink() {
  updateUrl(els.searchInput.value.trim());
  try {
    await navigator.clipboard.writeText(window.location.href);
    els.dataStatus.textContent = "Link copied";
  } catch {
    els.dataStatus.textContent = "Copy failed";
  }
}

function applyParams() {
  const params = new URLSearchParams(window.location.search);
  els.searchInput.value = params.get("q") || "";
  if (els.landingSearchInput) els.landingSearchInput.value = els.searchInput.value;
  for (const [param, el] of [
    ["trait", els.traitFilter],
    ["cohort", els.cohortFilter],
    ["ancestry", els.ancestryFilter],
    ["encoding", els.encodingFilter],
    ["p", els.pFilter],
  ]) {
    if (params.has(param)) el.value = params.get(param);
  }
  if (params.has("annotation")) {
    setCheckboxValues(els.annotationFilter, params.get("annotation").split(","));
  }
  const focusGene = params.get("focus_gene");
  const focusTrait = params.get("focus_trait");
  state.pendingFocus = focusGene && focusTrait ? { geneId: focusGene, trait: focusTrait } : null;
}

function attachEvents() {
  if (els.landingSearchForm) {
    els.landingSearchForm.addEventListener("submit", (event) => {
      event.preventDefault();
      els.searchInput.value = els.landingSearchInput.value.trim();
      runSearch();
    });
  }
  if (els.homeButton) {
    els.homeButton.addEventListener("click", () => {
      if (els.landingSearchInput) els.landingSearchInput.value = els.searchInput.value;
      showLandingPage();
      window.history.replaceState({}, "", window.location.pathname);
    });
  }
  els.searchButton.addEventListener("click", runSearch);
  els.topButton.addEventListener("click", renderTopHits);
  els.resetButton.addEventListener("click", resetFilters);
  els.downloadButton.addEventListener("click", downloadRows);
  if (els.downloadPlotButton) els.downloadPlotButton.addEventListener("click", downloadPlotPng);
  if (els.clearCohortFocusButton) els.clearCohortFocusButton.addEventListener("click", clearCohortFocus);
  els.copyLinkButton.addEventListener("click", copyLink);
  els.searchInput.addEventListener("keydown", (event) => {
    if (event.key === "Enter") runSearch();
  });
  for (const el of [
    els.traitFilter,
    els.cohortFilter,
    els.ancestryFilter,
    els.annotationFilter,
    els.encodingFilter,
    els.pFilter,
  ]) {
    el.addEventListener("change", () => {
      const base = state.lastLoadedRows.length ? state.lastLoadedRows : state.topHits;
      renderRows(base, els.resultTitle.textContent, "Filtered current rows.");
      updateUrl(els.searchInput.value.trim(), { focus: false });
    });
  }
}

async function init() {
  for (const [key, id] of Object.entries({
    landingPage: "landing-page",
    resultsPage: "results-page",
    landingEyebrow: "landing-eyebrow",
    landingTitle: "landing-title",
    landingCopy: "landing-copy",
    paperCitation: "paper-citation",
    paperPrefix: "paper-prefix",
    feedbackButton: "feedback-button",
    siteFootnote: "site-footnote",
    resultsFootnote: "results-footnote",
    landingSearchForm: "landing-search-form",
    landingSearchInput: "landing-search-input",
    searchSuggestions: "search-suggestions",
    homeButton: "home-button",
    dataStatus: "data-status",
    searchInput: "search-input",
    searchButton: "search-button",
    topButton: "top-button",
    traitFilter: "trait-filter",
    cohortFilter: "cohort-filter",
    ancestryFilter: "ancestry-filter",
    annotationFilter: "annotation-filter",
    encodingFilter: "encoding-filter",
    pFilter: "p-filter",
    resetButton: "reset-button",
    qcSummary: "qc-summary",
    summaryGrid: "summary-grid",
    resultTitle: "result-title",
    resultSubtitle: "result-subtitle",
    message: "message",
    metaTable: "meta-results-table",
    cohortTable: "cohort-results-table",
    cohortFocusLabel: "cohort-focus-label",
    clearCohortFocusButton: "clear-cohort-focus-button",
    detailPanel: "detail-panel",
    plotTitle: "plot-title",
    plotPanel: "plot-panel",
    downloadPlotButton: "download-plot-button",
    copyLinkButton: "copy-link-button",
    downloadButton: "download-button",
  })) {
    const element = $(id);
    els[key] = element;
  }

  try {
    [state.manifest, state.phenotypes, state.genes, state.symbolIndex, state.topHits, state.qc, state.siteContent] =
      await Promise.all([
        fetchJson(`${DATA_ROOT}/manifest.json`),
        fetchJson(`${DATA_ROOT}/phenotypes.json`),
        fetchJson(`${DATA_ROOT}/genes.json`),
        fetchJson(`${DATA_ROOT}/symbol_index.json`),
        fetchJson(`${DATA_ROOT}/top_hits.json`),
        fetchJson(`${DATA_ROOT}/qc/summary.json`).catch(() => []),
        fetchJson(`${DATA_ROOT}/site_content.json`).catch(() => null),
      ]);
  } catch (error) {
    els.dataStatus.textContent = "No generated data";
    showMessage(`Generated data files were not found. Run scripts/build_browser_data.py first. ${error.message}`);
    renderTable("meta", []);
    renderTable("cohort", []);
    return;
  }

  state.phenotypeById = new Map(state.phenotypes.map((item) => [item.id, item]));
  state.geneById = new Map(state.genes.map((item) => [item.id, item]));
  setupFilters();
  populateSearchSuggestions();
  applySiteContent();
  renderSummary();
  renderQc();
  attachEvents();
  applyParams();
  state.lastLoadedRows = state.topHits;
  els.dataStatus.textContent = `${state.topHits.length.toLocaleString()} indexed hits`;
  const hasResultParams = new URLSearchParams(window.location.search).toString() !== "";
  if (els.searchInput.value.trim()) {
    await runSearch({ preserveFocus: true });
  } else if (hasResultParams) {
    showResultsPage();
    renderRows(state.topHits, "Top associations", "Smallest available p-values from the generated index.");
  } else {
    showLandingPage();
  }
}

init();

