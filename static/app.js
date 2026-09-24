const $ = id => document.getElementById(id);
let offset = 0, gene = '', selectedId = null;
const LIMIT = 25;
const value = x => x === null || x === undefined || x === '' ? 'Not available' : String(x);
const frequency = n => n === null || n === undefined ? 'Not available' : Number(n).toExponential(2);

function add(parent, tag, text, className) {
  const el = document.createElement(tag);
  el.textContent = text;
  if (className) el.className = className;
  parent.append(el);
  return el;
}
function field(parent, key, val) {
  const row = add(parent, 'div', '', 'field');
  add(row, 'span', key); add(row, 'span', value(val));
}

async function list() {
  $('status').textContent = 'Loading…';
  try {
    const q = new URLSearchParams({limit: LIMIT, offset, gene});
    const response = await fetch(`/api/variants?${q}`);
    if (!response.ok) throw Error(`API returned ${response.status}`);
    const {items} = await response.json();
    const body = $('variants'); body.replaceChildren();
    for (const item of items) {
      const tr = add(body, 'tr', '');
      tr.tabIndex = 0;
      tr.setAttribute('aria-label', `${item.chrom}:${item.pos} ${item.ref} to ${item.alt}`);
      if (item.id === selectedId) tr.className = 'selected';
      for (const cell of [`${item.chrom}:${item.pos}`, `${item.ref} → ${item.alt}`,
        value(item.gene_symbol), value(item.consequence).replaceAll('_', ' '), frequency(item.gnomadg_af)]) add(tr, 'td', cell);
      tr.addEventListener('click', () => show(item.id));
      tr.addEventListener('keydown', e => { if (e.key === 'Enter') show(item.id); });
    }
    $('status').textContent = items.length ? '' : 'No variants found. Run an annotation import or change the gene filter.';
    $('count').textContent = `${items.length} on this page`;
    $('page').textContent = `Page ${Math.floor(offset / LIMIT) + 1}`;
    $('prev').disabled = offset === 0;
    $('next').disabled = items.length < LIMIT;
  } catch (error) { $('status').textContent = `Could not load variants: ${error.message}`; }
}

async function show(id) {
  selectedId = id;
  document.querySelectorAll('tbody tr').forEach(tr => tr.classList.remove('selected'));
  const panel = $('detail'); panel.replaceChildren(); add(panel, 'h2', 'Loading detail…');
  try {
    const response = await fetch(`/api/variants/${id}`);
    if (!response.ok) throw Error(`API returned ${response.status}`);
    const v = await response.json();
    panel.replaceChildren();
    add(panel, 'h2', `${v.chrom}:${v.pos} ${v.ref} → ${v.alt}`);
    add(panel, 'span', 'GRCh38 · VEP annotation', 'badge');
    field(panel, 'Gene', v.gene_symbol);
    field(panel, 'Gene ID', v.gene_id);
    field(panel, 'Selected transcript', v.transcript_id);
    field(panel, 'Most severe consequence', v.consequence);
    field(panel, 'gnomAD genomes AF', frequency(v.gnomadg_af));
    add(panel, 'h3', 'Data provenance');
    field(panel, 'Source', v.source);
    field(panel, 'Version', v.source_version);
    field(panel, 'Imported', new Date(v.annotated_at).toLocaleString());
    add(panel, 'h3', 'ClinVar condition record');
    add(panel, 'p', 'This is a dated assertion from ClinVar, imported separately from VEP. It is not a classification made by this application.', 'note');
    if (v.clinvar_assertions.length) {
      for (const a of v.clinvar_assertions) {
        field(panel, 'Condition', `${a.condition} (${a.disease_id})`);
        field(panel, 'Reported classification', a.classification);
        field(panel, 'Review status', a.review_status);
        field(panel, 'Record snapshot', `${a.record_accession} · ${a.snapshot_date}`);
        const link = add(panel, 'a', 'View condition record on ClinVar');
        link.href = a.source_url;
        link.target = '_blank';
        link.rel = 'noopener noreferrer';
      }
    } else add(panel, 'p', 'No ClinVar condition record imported for this allele.');
    add(panel, 'h3', 'Gene–phenotype associations');
    add(panel, 'p', 'These HPO associations belong to a gene and disease, not to this variant. They cannot establish causality.', 'note');
    if (v.hpo_condition_ids.length) add(panel, 'p', `Showing HPO associations for ${v.hpo_condition_ids.join(', ')}.`);
    if (v.gene_phenotypes.length) {
      const ul = add(panel, 'ul', '', 'hpo-list');
      for (const a of v.gene_phenotypes) add(ul, 'li', `${a.hpo_id}: ${a.name} (${a.disease_id})`);
      field(panel, 'HPO source SHA-256', v.gene_phenotypes[0].source_sha256);
      if (v.gene_phenotypes.length === 250) add(panel, 'p', 'Showing first 250 associations.');
    } else add(panel, 'p', 'No imported HPO associations for the selected gene.');
    add(panel, 'h3', 'Full annotation');
    const details = add(panel, 'details', '');
    add(details, 'summary', 'Show raw VEP JSON');
    add(details, 'pre', JSON.stringify(v.vep_json, null, 2));
  } catch (error) { panel.replaceChildren(); add(panel, 'p', `Could not load detail: ${error.message}`); }
}

$('search').addEventListener('submit', e => { e.preventDefault(); gene = $('gene').value.trim(); offset = 0; list(); });
$('prev').addEventListener('click', () => { offset = Math.max(0, offset - LIMIT); list(); });
$('next').addEventListener('click', () => { offset += LIMIT; list(); });
list();
