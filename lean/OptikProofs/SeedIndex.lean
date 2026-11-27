import OptikProofs.Minimizer

/-!
# Parallel seed extraction equals the serial scan (prop:seed)

This is step 1 of `SeedIndex::build` in `include/seed_index.hpp`.
Each task is assigned the positions [s, e) and scans the windows that start in [s - (w - 1), e), clipped at 0.
The last task of a reference takes every position from s on.
A task reports a window's minimizer only when the position is its own and differs from the last one reported, as `extract_minimizers` does with `last_pos`.

The theorem says the task outputs, concatenated in task order, are exactly the serial scan's output.
So every seed position appears once and the postings are in reference order, whatever the cut points.
The proof rests on two facts.
Minimizer positions never decrease from one window to the next (`Minimizer.minPos_mono`), so both outputs are strictly increasing lists.
A window can only select a position within w - 1 of its start, so a task's windows include every window that selects one of its positions.
-/

namespace OptikProofs.SeedIndex

open OptikProofs.Minimizer

/-- `extract_minimizers` emits a position unless it equals the last one emitted. -/
def emitNew : Option ℕ → List ℕ → List ℕ
  | _, [] => []
  | last, p :: ps => if last = some p then emitNew last ps else p :: emitNew (some p) ps

/-- On a non-decreasing input whose values are all at least `x`, emitting after `x` yields a strictly increasing list of the input's values other than `x`. -/
theorem emitNew_some_spec :
    ∀ (l : List ℕ) (x : ℕ), l.Pairwise (· ≤ ·) → (∀ y ∈ l, x ≤ y) →
      (∀ a, a ∈ emitNew (some x) l ↔ a ∈ l ∧ a ≠ x) ∧ (emitNew (some x) l).Pairwise (· < ·) ∧
        ∀ y ∈ emitNew (some x) l, x < y
  | [], x, _, _ => by simp [emitNew]
  | p :: ps, x, hs, hx => by
    rw [List.pairwise_cons] at hs
    by_cases hxp : x = p
    · subst hxp
      obtain ⟨hmem, hpw, hgt⟩ := emitNew_some_spec ps x hs.2 hs.1
      rw [emitNew, if_pos rfl]
      refine ⟨fun a => ?_, hpw, hgt⟩
      rw [hmem]
      simp only [List.mem_cons]
      constructor
      · rintro ⟨h1, h2⟩; exact ⟨Or.inr h1, h2⟩
      · rintro ⟨h1 | h1, h2⟩
        · exact absurd h1 h2
        · exact ⟨h1, h2⟩
    · have hxlt : x < p := lt_of_le_of_ne (hx p (by simp)) hxp
      obtain ⟨hmem, hpw, hgt⟩ := emitNew_some_spec ps p hs.2 hs.1
      simp only [emitNew, Option.some.injEq, if_neg hxp]
      refine ⟨fun a => ?_, List.pairwise_cons.mpr ⟨hgt, hpw⟩, fun y hy => ?_⟩
      · simp only [List.mem_cons, hmem]
        constructor
        · rintro (rfl | ⟨h1, -⟩)
          · exact ⟨Or.inl rfl, fun h => hxp h.symm⟩
          · exact ⟨Or.inr h1, fun h => by subst h; exact absurd (hs.1 a h1) (not_le.mpr hxlt)⟩
        · rintro ⟨rfl | h1, h2⟩
          · exact Or.inl rfl
          · by_cases hap : a = p
            · exact Or.inl hap
            · exact Or.inr ⟨h1, hap⟩
      · simp only [List.mem_cons] at hy
        rcases hy with rfl | hy
        · exact hxlt
        · exact lt_trans hxlt (hgt y hy)

/-- From the start, with nothing emitted yet, the output is the input's values in strictly increasing order. -/
theorem emitNew_none_spec (l : List ℕ) (hs : l.Pairwise (· ≤ ·)) :
    (∀ a, a ∈ emitNew none l ↔ a ∈ l) ∧ (emitNew none l).Pairwise (· < ·) := by
  cases l with
  | nil => simp [emitNew]
  | cons p ps =>
    rw [List.pairwise_cons] at hs
    obtain ⟨hmem, hpw, hgt⟩ := emitNew_some_spec ps p hs.2 hs.1
    simp only [emitNew, reduceCtorEq, if_false]
    refine ⟨fun a => ?_, List.pairwise_cons.mpr ⟨hgt, hpw⟩⟩
    simp only [List.mem_cons, hmem]
    constructor
    · rintro (rfl | ⟨h, -⟩)
      · exact Or.inl rfl
      · exact Or.inr h
    · rintro (rfl | h)
      · exact Or.inl rfl
      · by_cases hap : a = p
        · exact Or.inl hap
        · exact Or.inr ⟨h, hap⟩

/-- Two strictly increasing lists with the same members are the same list. -/
theorem eq_of_strictMono_mem {l₁ l₂ : List ℕ} (h₁ : l₁.Pairwise (· < ·)) (h₂ : l₂.Pairwise (· < ·))
    (hmem : ∀ a, a ∈ l₁ ↔ a ∈ l₂) : l₁ = l₂ := by
  have hperm : l₁.Perm l₂ :=
    (List.perm_ext_iff_of_nodup (h₁.imp ne_of_lt) (h₂.imp ne_of_lt)).mpr hmem
  exact List.Perm.eq_of_pairwise' (r := (· ≤ ·)) (h₁.imp le_of_lt) (h₂.imp le_of_lt) hperm

variable {α : Type*} [LinearOrder α] (h : ℕ → α) (w nwin : ℕ)

/-- Window minimizer positions over a list of window starts. -/
def mins (starts : List ℕ) : List ℕ := starts.map (minPos h w)

/-- The serial scan over all `nwin` windows of the reference. -/
def serial : List ℕ := emitNew none (mins h w (List.range nwin))

/-- The windows a task assigned [s, e) scans, starting at s - (w - 1), which natural subtraction clips at 0. -/
def taskWindows (s e : ℕ) : List ℕ := List.range' (s - (w - 1)) (e - (s - (w - 1)))

/-- Whether position p belongs to a task starting at s, bounded by `ownEnd` (the last task has none). -/
def owns (s : ℕ) (ownEnd : Option ℕ) (p : ℕ) : Bool :=
  decide (s ≤ p) && match ownEnd with
    | none => true
    | some e => decide (p < e)

/-- One task's output, which filters its windows' minimizers to its own positions and reports them as `extract_minimizers` does. -/
def taskOut (s e : ℕ) (ownEnd : Option ℕ) : List ℕ :=
  emitNew none ((mins h w (taskWindows w s e)).filter (owns s ownEnd))

/-- Tasks from window start `s` through the cut points `cuts`, concatenated in order.
Each task ends at the next cut and takes the positions up to it, while the last runs to `nwin` and takes everything after s. -/
def runTasks : ℕ → List ℕ → List ℕ
  | s, [] => taskOut h w s nwin none
  | s, c :: cs => taskOut h w s c (some c) ++ runTasks c cs

variable {h w nwin}

theorem mins_sorted (hw : 1 ≤ w) {starts : List ℕ} (hs : starts.Pairwise (· < ·)) :
    (mins h w starts).Pairwise (· ≤ ·) :=
  List.pairwise_map.mpr (hs.imp fun hab => minPos_monotone h hw (le_of_lt hab))

theorem mem_serial (hw : 1 ≤ w) (a : ℕ) : a ∈ serial h w nwin ↔ ∃ j < nwin, minPos h w j = a := by
  rw [serial, (emitNew_none_spec _ (mins_sorted hw List.pairwise_lt_range)).1, mins, List.mem_map]
  simp

theorem serial_sorted (hw : 1 ≤ w) : (serial h w nwin).Pairwise (· < ·) :=
  (emitNew_none_spec _ (mins_sorted hw List.pairwise_lt_range)).2

/-- A task finds every minimizer position assigned to it.
A window selects a position at most w - 1 past its start, so a window choosing one of the task's positions starts at s - (w - 1) or later, and before the task's end. -/
theorem mem_taskOut (hw : 1 ≤ w) {s e : ℕ} {ownEnd : Option ℕ} (he : e ≤ nwin)
    (hown : ownEnd = some e ∨ (ownEnd = none ∧ e = nwin)) (a : ℕ) :
    a ∈ taskOut h w s e ownEnd ↔ (∃ j < nwin, minPos h w j = a) ∧ owns s ownEnd a = true := by
  have hsorted : ((mins h w (taskWindows w s e)).filter (owns s ownEnd)).Pairwise (· ≤ ·) :=
    (mins_sorted hw List.pairwise_lt_range').filter _
  rw [taskOut, (emitNew_none_spec _ hsorted).1, List.mem_filter, mins, List.mem_map]
  simp only [taskWindows, List.mem_range']
  constructor
  · rintro ⟨⟨j, ⟨i, hi, rfl⟩, rfl⟩, ho⟩
    exact ⟨⟨_, by omega, rfl⟩, ho⟩
  · rintro ⟨⟨j, hj, rfl⟩, ho⟩
    obtain ⟨hlo, hhi, -, -⟩ := minPos_spec h hw j
    have hs : s ≤ minPos h w j := by
      simp only [owns, Bool.and_eq_true, decide_eq_true_eq] at ho
      exact ho.1
    have hje : j < e := by
      rcases hown with rfl | ⟨rfl, rfl⟩
      · simp only [owns, Bool.and_eq_true, decide_eq_true_eq] at ho
        omega
      · exact hj
    exact ⟨⟨j, ⟨j - (s - (w - 1)), by omega, by omega⟩, rfl⟩, ho⟩

theorem taskOut_sorted (hw : 1 ≤ w) (s e : ℕ) (ownEnd : Option ℕ) :
    (taskOut h w s e ownEnd).Pairwise (· < ·) :=
  (emitNew_none_spec _ ((mins_sorted hw List.pairwise_lt_range').filter _)).2

/-- Tasks from `s` through the cuts report, in increasing order, every minimizer position at or after `s`. -/
theorem runTasks_spec (hw : 1 ≤ w) :
    ∀ (cuts : List ℕ) (s : ℕ), (s :: cuts).Pairwise (· < ·) → (∀ c ∈ cuts, c ≤ nwin) → s ≤ nwin →
      (∀ a, a ∈ runTasks h w nwin s cuts ↔ (∃ j < nwin, minPos h w j = a) ∧ s ≤ a) ∧
        (runTasks h w nwin s cuts).Pairwise (· < ·)
  | [], s, _, _, _ => by
    refine ⟨fun a => ?_, taskOut_sorted hw _ _ _⟩
    rw [runTasks, mem_taskOut hw le_rfl (Or.inr ⟨rfl, rfl⟩)]
    simp [owns]
  | c :: cs, s, hsorted, hcuts, _ => by
    rw [List.pairwise_cons] at hsorted
    have hsc : s < c := hsorted.1 c (by simp)
    have hc : c ≤ nwin := hcuts c (by simp)
    obtain ⟨hmem, hpw⟩ := runTasks_spec hw cs c hsorted.2 (fun c' hc' => hcuts c' (by simp [hc'])) hc
    have hfirst := mem_taskOut (h := h) hw (s := s) hc (Or.inl rfl)
    refine ⟨fun a => ?_, ?_⟩
    · rw [runTasks, List.mem_append, hfirst, hmem]
      simp only [owns, Bool.and_eq_true, decide_eq_true_eq]
      constructor
      · rintro (⟨hj, hs, -⟩ | ⟨hj, hc'⟩)
        · exact ⟨hj, hs⟩
        · exact ⟨hj, by omega⟩
      · rintro ⟨hj, hs⟩
        by_cases hac : a < c
        · exact Or.inl ⟨hj, hs, hac⟩
        · exact Or.inr ⟨hj, by omega⟩
    · rw [runTasks, List.pairwise_append]
      refine ⟨taskOut_sorted hw _ _ _, hpw, fun x hx y hy => ?_⟩
      have hx' := ((hfirst x).mp hx).2
      simp only [owns, Bool.and_eq_true, decide_eq_true_eq] at hx'
      have hy' := ((hmem y).mp hy).2
      omega

/-- prop:seed.
For any strictly increasing cut points inside the reference, the concatenated task outputs equal the serial scan.
Every seed position therefore appears exactly once, with postings in reference order. -/
theorem parallel_eq_serial (hw : 1 ≤ w) (cuts : List ℕ) (hsorted : (0 :: cuts).Pairwise (· < ·))
    (hcuts : ∀ c ∈ cuts, c ≤ nwin) : runTasks h w nwin 0 cuts = serial h w nwin := by
  obtain ⟨hmem, hpw⟩ := runTasks_spec (h := h) hw cuts 0 hsorted hcuts (Nat.zero_le _)
  refine eq_of_strictMono_mem hpw (serial_sorted (h := h) hw) fun a => ?_
  rw [hmem, mem_serial (h := h) hw]
  simp

end OptikProofs.SeedIndex
