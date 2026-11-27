import Mathlib.Data.Finset.Basic

/-!
# An editable k-mer set over a static core

Models `DynSpss` in `include/dyn_spss.hpp`, which `spss-update` drives through `okdyn_apply_fasta` and `okdyn_compact_save` in `src/ffi_spss.cpp`.

Keys are canonical k-mers.
`core` is the static index's set, whose `contains` answers exactly membership of the canonical form by `Spss.present_iff`.
`added` and `deleted` are the two DeltaOA sets, modelled as finite sets because DeltaOA's own behaviour is a hash-table question outside this model.
The set the structure stands for is `(core \ deleted) ∪ added`.

* After `insert` or `erase`, `added` is still disjoint from the core and `deleted` still inside it.
* Each one changes the set by exactly that one k-mer.
* Their return values say whether membership flipped, which is what `okdyn_apply_fasta` adds up as "changed".
* `compact` builds its new core from exactly the current set, so saving loses and invents nothing.
-/

namespace OptikProofs.Dyn

/-- The core's k-mers and the two edit buffers. -/
structure State where
  core : Finset ℕ
  added : Finset ℕ
  deleted : Finset ℕ

/-- `added` holds only k-mers the core lacks, and `deleted` only k-mers it has. -/
def State.Wf (st : State) : Prop := Disjoint st.added st.core ∧ st.deleted ⊆ st.core

/-- The k-mer set the structure stands for. -/
def State.set (st : State) : Finset ℕ := (st.core \ st.deleted) ∪ st.added

/-- `DynSpss::insert`, returning the new state and whether the k-mer was absent before. -/
def insert (st : State) (c : ℕ) : State × Bool :=
  if c ∈ st.core then ({ st with deleted := st.deleted.erase c }, decide (c ∈ st.deleted))
  else ({ st with added := Insert.insert c st.added }, decide (c ∉ st.added))

/-- `DynSpss::erase`, which also reports whether the k-mer was present before. -/
def erase (st : State) (c : ℕ) : State × Bool :=
  if c ∈ st.core then ({ st with deleted := Insert.insert c st.deleted }, decide (c ∉ st.deleted))
  else ({ st with added := st.added.erase c }, decide (c ∈ st.added))

/-- The k-mer list `compact` hands to the graph build, before sorting and deduplication. -/
def compactInput (st : State) : Finset ℕ := st.core.filter (· ∉ st.deleted) ∪ st.added

/-- After a successful `compact`, the new core holds the deduplicated input, and both buffers are empty. -/
def compact (st : State) : State := ⟨compactInput st, ∅, ∅⟩

theorem insert_wf {st : State} (h : st.Wf) (c : ℕ) : (insert st c).1.Wf := by
  unfold insert
  split_ifs with hc
  · exact ⟨h.1, (Finset.erase_subset _ _).trans h.2⟩
  · refine ⟨?_, h.2⟩
    exact Finset.disjoint_insert_left.mpr ⟨hc, h.1⟩

theorem erase_wf {st : State} (h : st.Wf) (c : ℕ) : (erase st c).1.Wf := by
  unfold erase
  split_ifs with hc
  · exact ⟨h.1, Finset.insert_subset hc h.2⟩
  · exact ⟨Finset.disjoint_of_subset_left (Finset.erase_subset _ _) h.1, h.2⟩

theorem mem_set {st : State} (c : ℕ) : c ∈ st.set ↔ (c ∈ st.core ∧ c ∉ st.deleted) ∨ c ∈ st.added := by
  simp [State.set]

/-- Insert adds exactly that k-mer. -/
theorem insert_set (st : State) (c : ℕ) : (insert st c).1.set = Insert.insert c st.set := by
  ext x
  unfold insert
  split_ifs with hc
  · simp only [mem_set, Finset.mem_insert, Finset.mem_erase]
    by_cases hx : x = c
    · subst hx
      simp [hc]
    · simp [hx]
  · simp only [mem_set, Finset.mem_insert]
    by_cases hx : x = c
    · subst hx
      simp
    · simp [hx]

/-- Erase removes exactly that k-mer. -/
theorem erase_set {st : State} (h : st.Wf) (c : ℕ) : (erase st c).1.set = st.set.erase c := by
  ext x
  unfold erase
  split_ifs with hc
  · simp only [mem_set, Finset.mem_insert, Finset.mem_erase]
    by_cases hx : x = c
    · subst hx
      have : x ∉ st.added := fun ha => Finset.disjoint_left.mp h.1 ha hc
      simp [this]
    · simp [hx]
  · simp only [mem_set, Finset.mem_erase]
    by_cases hx : x = c
    · subst hx
      simp [hc]
    · simp [hx]

/-- Insert reports true exactly when the k-mer was missing. -/
theorem insert_result {st : State} (h : st.Wf) (c : ℕ) : (insert st c).2 = decide (c ∉ st.set) := by
  unfold insert
  split_ifs with hc
  · have : c ∉ st.added := fun ha => Finset.disjoint_left.mp h.1 ha hc
    simp [mem_set, hc, this]
  · have : c ∉ st.deleted := fun hd => hc (h.2 hd)
    simp [mem_set, hc]

/-- Erase reports true exactly when the k-mer was there. -/
theorem erase_result {st : State} (h : st.Wf) (c : ℕ) : (erase st c).2 = decide (c ∈ st.set) := by
  unfold erase
  split_ifs with hc
  · have : c ∉ st.added := fun ha => Finset.disjoint_left.mp h.1 ha hc
    simp [mem_set, hc, this]
  · simp [mem_set, hc]

/-- `compact` rebuilds exactly the current set, and the new state again satisfies the invariants. -/
theorem compact_set (st : State) : (compact st).set = st.set ∧ (compact st).Wf := by
  refine ⟨?_, by simp [compact, State.Wf]⟩
  ext x
  simp [compact, compactInput, State.set]

end OptikProofs.Dyn
