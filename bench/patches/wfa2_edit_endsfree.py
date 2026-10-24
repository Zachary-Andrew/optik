#!/usr/bin/env python3
"""wfa2_edit_endsfree.py - idempotent patch that makes ends-free (semi-global) alignment work under the edit and indel metrics of the vendored WFA2-lib.

Why this patch exists
optik's mapper contract scores reads semi-globally (query fully consumed, reference ends free) with unit-cost edit distance (match=0, mismatch=1, indel=1).
WFA2-lib computes ends-free alignment with the "boundary-key" formula `score / -penalties->match`, and bails out of ends-free initialization whenever `penalties->match == 0`.
The edit and indel metrics set match=0, so ends-free alignment was silently a no-op for the exact metric the mapper needs, and every "ends-free" edit alignment came back as plain end-to-end.

The boundary key is only an implementation artifact.
For unit-cost metrics every edit operation, each freely skipped base included, advances the score by exactly 1, so the correct boundary key is `score` itself.
For the linear and affine family it stays `score / -match`.
This patch substitutes that key at the four ends-free sites in wavefront_compute.c and relaxes the match==0 guard for edit/indel only.

build_tools.sh build_wfa2() applies it AFTER clone_or_update hard-resets to origin/HEAD, so it always runs against a pristine tree.
It is idempotent.
If the marker function already exists the file is left untouched, and if upstream drift breaks any anchor line it exits 1 without modifying the file, so the build can warn and continue.

Usage:  wfa2_edit_endsfree.py <path/to/wavefront_compute.c> Exit:   0 when patched or already patched, 1 on an anchor mismatch with the file untouched
"""

import sys

MARKER = "wavefront_compute_endsfree_key"

# One replacement of the whole endsfree_required preamble, which prepends the boundary-key helper and rewrites the guard and divisor logic.
REQUIRED_BEFORE = """bool wavefront_compute_endsfree_required(
    wavefront_aligner_t* const wf_aligner,
    const int score) {
  // Parameters
  alignment_form_t* const alg_form = &wf_aligner->alignment_form;
  wavefront_penalties_t* const penalties = &wf_aligner->penalties;
  // Return if ends-free initialization is required
  if (penalties->match == 0) return false;
  if (alg_form->span != alignment_endsfree) return false;
  if (alg_form->text_begin_free == 0 &&
      alg_form->pattern_begin_free == 0) return false;
  if (score % (-penalties->match) != 0) return false;
  // Check boundary conditions for ends-free
  const int endsfree_k = score/(-penalties->match); // (h/v)-coordinate for boundary conditions
  const bool text_begin_free = (alg_form->text_begin_free >= endsfree_k);
  const bool pattern_begin_free = (alg_form->pattern_begin_free >= endsfree_k);
  if (!text_begin_free && !pattern_begin_free) return false;
  // Ok
  return true;
}"""

REQUIRED_AFTER = """/*
 * Ends-free boundary key
 *
 * Maps a score s to the diagonal of the ends-free boundary point
 * ((h=s/div, v=0) for text-begin-free, (h=0, v=s/div) for pattern-begin-free).
 *
 * The linear/affine family uses the transformed (Eizenga) units where the
 * boundary advances one diagonal per -match score units, so div = -match.
 * The edit/indel metrics are unit-cost (match == 0, and the original code
 * bailed on match == 0 because it divided by -match): every edit operation
 * -- including every free-skipped base -- moves the score by exactly 1, so
 * the boundary advances one diagonal per score unit, div = 1.  This is the
 * fix that lets ends-free (semi-global) alignment work under the edit
 * distance metric (optik's mapper contract).
 */
static inline int wavefront_compute_endsfree_key(
    const wavefront_penalties_t* const penalties,
    const int score) {
  const bool unit_cost = (penalties->match == 0 &&
      (penalties->distance_metric == edit ||
       penalties->distance_metric == indel));
  const int div = unit_cost ? 1 : -penalties->match; // caller guards div>0
  return score / div;
}
bool wavefront_compute_endsfree_required(
    wavefront_aligner_t* const wf_aligner,
    const int score) {
  // Parameters
  alignment_form_t* const alg_form = &wf_aligner->alignment_form;
  wavefront_penalties_t* const penalties = &wf_aligner->penalties;
  // Return if ends-free initialization is required
  if (penalties->match == 0) {
    // Unit-cost metrics (edit/indel) are the only match==0 metrics whose
    // score unit is 1 (linear/affine with match==0 advance by their own
    // mismatch/indel units and have no defined boundary here).
    if (penalties->distance_metric != edit &&
        penalties->distance_metric != indel) return false;
  } else {
    if (score % (-penalties->match) != 0) return false;
  }
  if (alg_form->span != alignment_endsfree) return false;
  if (alg_form->text_begin_free == 0 &&
      alg_form->pattern_begin_free == 0) return false;
  // Check boundary conditions for ends-free
  const int endsfree_k = wavefront_compute_endsfree_key(penalties,score); // (h/v)-coordinate for boundary conditions
  const bool text_begin_free = (alg_form->text_begin_free >= endsfree_k);
  const bool pattern_begin_free = (alg_form->pattern_begin_free >= endsfree_k);
  if (!text_begin_free && !pattern_begin_free) return false;
  // Ok
  return true;
}"""

# The remaining three sites use the same divisor expression.
CONST_KEY_BEFORE = "  const int endsfree_k = score/(-penalties->match);"
CONST_KEY_AFTER = (
    "  const int endsfree_k = wavefront_compute_endsfree_key(penalties,score);"
)
INT_KEY_BEFORE = "  int endsfree_k = score/(-penalties->match);"
INT_KEY_AFTER = "  int endsfree_k = wavefront_compute_endsfree_key(penalties,score);"


def main() -> int:
    """Patch the file named on the command line, returning the exit status."""
    if len(sys.argv) != 2:
        print("usage: " + __doc__.split("Usage:")[1].strip(), file=sys.stderr)
        return 2
    path = sys.argv[1]
    with open(path, "r", encoding="utf-8") as f:
        text = f.read()

    if MARKER in text:
        print(f"{path}: already patched (marker {MARKER!r} present), no-op")
        return 0

    # The required() anchor must be present exactly once.
    if text.count(REQUIRED_BEFORE) != 1:
        print(
            f"{path}: anchor for endsfree_required not found exactly once "
            f"(count={text.count(REQUIRED_BEFORE)}); upstream drift? no changes made",
            file=sys.stderr,
        )
        return 1

    new_text = text.replace(REQUIRED_BEFORE, REQUIRED_AFTER)
    # The const-key expression appears twice, in endsfree_limits and endsfree_allocate_null, and the int-key once, in endsfree_init.
    ck = new_text.count(CONST_KEY_BEFORE)
    ik = new_text.count(INT_KEY_BEFORE)
    if ck != 2 or ik != 1:
        print(
            f"{path}: key-expression anchors shifted (const={ck} int={ik}); "
            "upstream drift? no changes made",
            file=sys.stderr,
        )
        return 1

    new_text = new_text.replace(CONST_KEY_BEFORE, CONST_KEY_AFTER)
    new_text = new_text.replace(INT_KEY_BEFORE, INT_KEY_AFTER)

    if MARKER not in new_text:
        print(
            f"{path}: patch produced no marker; aborting (no changes made)",
            file=sys.stderr,
        )
        return 1

    with open(path, "w", encoding="utf-8") as f:
        f.write(new_text)
    print(f"{path}: patched (edit/indel ends-free boundary key)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
