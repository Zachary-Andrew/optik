import OptikProofs.Align

/-!
# CIGAR bookkeeping

Models `push_op`, `ref_span`, `cigar_cols`, `cigar_ed_implied`, `clip_ends` and `hard_clip` in `tools/src/map/cigar.rs` and `tools/src/map/paired.rs`.

A CIGAR is a list of runs, each a length and an op.
Parsing the digits is left to the test vectors, which check the code against `edImplied`, `refSpan` and `colCount` here.

* `push_op` appends exactly the columns asked for, and never leaves an empty run or two neighbouring runs with the same op (`expand_pushOp`, `canon_pushOp`).
* `ref_span`, which places END and TLEN, counts exactly the reference bases the columns consume, and `cigar_cols` counts the columns (`refSpan_eq`, `colCount_eq`).
* The NM that `cigar_ed_implied` recomputes is `Align.cost` of the CIGAR's columns, the cost `Align.shift_deletion` and `Align.shift_insertion` are about (`edImplied_eq_cost`).
* `clip_ends` preserves the read's length, its new offset is where the core starts, and its NM is the core's edit cost, never more than the whole alignment's (`clip_query_length`, `clip_nm`, `clip_nm_le`).
  After `hard_clip` the CIGAR's query length equals the cut sequence's and the reference span is unchanged (`hardClip_qlen`, `hardClip_refSpan`).
-/

namespace OptikProofs.Cigar

/-- The SAM CIGAR operations optik writes or reads. -/
inductive COp
  | M | I | D | N | S | H | Eq | X
  deriving DecidableEq, Repr

/-- The ops that consume query bases, which SEQ must match. -/
def COp.query : COp → Bool
  | .M | .I | .S | .Eq | .X => true
  | _ => false

/-- The ops `ref_span` adds up. -/
def COp.ref : COp → Bool
  | .M | .D | .N | .Eq | .X => true
  | _ => false

/-- The ops `cigar_cols` adds up. -/
def COp.col : COp → Bool
  | .M | .Eq | .X | .I | .D => true
  | _ => false

/-- Run-length CIGAR. -/
abbrev Rle := List (ℕ × COp)

/-- One column per base of every run. -/
def expand (c : Rle) : List COp := c.flatMap fun r => List.replicate r.1 r.2

/-- The sum of the run lengths whose op satisfies `p`, as the code's digit loops add them up. -/
def runSum (p : COp → Bool) : Rle → ℕ
  | [] => 0
  | (n, o) :: rest => (if p o then n else 0) + runSum p rest

def refSpan : Rle → ℕ := runSum COp.ref
def colCount : Rle → ℕ := runSum COp.col
def qlen : Rle → ℕ := runSum COp.query

theorem runSum_eq (p : COp → Bool) : ∀ c : Rle, runSum p c = (expand c).countP p
  | [] => rfl
  | (n, o) :: rest => by
    simp only [runSum, expand, List.flatMap_cons, List.countP_append, List.countP_replicate] at *
    rw [runSum_eq p rest, expand]

theorem refSpan_eq (c : Rle) : refSpan c = (expand c).countP COp.ref := runSum_eq _ c
theorem colCount_eq (c : Rle) : colCount c = (expand c).countP COp.col := runSum_eq _ c
theorem qlen_eq (c : Rle) : qlen c = (expand c).countP COp.query := runSum_eq _ c

theorem runSum_append (p : COp → Bool) (a b : Rle) : runSum p (a ++ b) = runSum p a + runSum p b := by
  induction a with
  | nil => simp [runSum]
  | cons r a ih => obtain ⟨n, o⟩ := r; simp only [List.cons_append, runSum, ih]; omega

/-! ## push_op -/

/-- `push_op` appends `n` columns of `o`, merging with the last run when it has the same op. -/
def pushOp (n : ℕ) (o : COp) : Rle → Rle
  | [] => if n = 0 then [] else [(n, o)]
  | [(m, o')] => if n = 0 then [(m, o')] else if o' = o then [(m + n, o)] else [(m, o'), (n, o)]
  | x :: y :: r => x :: pushOp n o (y :: r)

/-- No empty run, and neighbouring runs differ in op. -/
def Canon : Rle → Prop
  | [] => True
  | [x] => 0 < x.1
  | x :: y :: r => 0 < x.1 ∧ x.2 ≠ y.2 ∧ Canon (y :: r)

theorem expand_pushOp (n : ℕ) (o : COp) : ∀ c : Rle, expand (pushOp n o c) = expand c ++ List.replicate n o
  | [] => by by_cases h : n = 0 <;> simp [pushOp, expand, h]
  | [(m, o')] => by
    by_cases h : n = 0
    · simp [pushOp, expand, h]
    · by_cases ho : o' = o
      · subst ho; simp [pushOp, expand, h]
      · simp [pushOp, expand, h, ho]
  | x :: y :: r => by
    have ih := expand_pushOp n o (y :: r)
    simp only [pushOp, expand, List.flatMap_cons] at ih ⊢
    rw [ih]
    simp

/-- The first run's op is unchanged. -/
theorem pushOp_head (n : ℕ) (o : COp) (y : ℕ × COp) (r : Rle) :
    ∃ z rest, pushOp n o (y :: r) = z :: rest ∧ z.2 = y.2 := by
  match r with
  | [] =>
    obtain ⟨m, o'⟩ := y
    by_cases h : n = 0
    · exact ⟨(m, o'), [], by simp [pushOp, h], rfl⟩
    · by_cases ho : o' = o
      · exact ⟨(m + n, o), [], by simp [pushOp, h, ho], ho.symm⟩
      · exact ⟨(m, o'), [(n, o)], by simp [pushOp, h, ho], rfl⟩
  | y' :: r' => exact ⟨y, _, rfl, rfl⟩

theorem canon_pushOp (n : ℕ) (o : COp) : ∀ c : Rle, Canon c → Canon (pushOp n o c)
  | [], _ => by by_cases h : n = 0 <;> simp [pushOp, Canon, h]; omega
  | [(m, o')], hc => by
    simp only [Canon] at hc
    by_cases h : n = 0
    · simp [pushOp, Canon, h, hc]
    · by_cases ho : o' = o
      · simp [pushOp, Canon, h, ho]; omega
      · simp [pushOp, Canon, h, ho, hc]; omega
  | x :: y :: r, ⟨hx, hxy, hc⟩ => by
    have ih := canon_pushOp n o (y :: r) hc
    obtain ⟨z, rest, he, hz⟩ := pushOp_head n o y r
    simp only [pushOp]
    rw [he] at ih ⊢
    cases rest with
    | nil => exact ⟨hx, hz ▸ hxy, ih⟩
    | cons w rest => exact ⟨hx, hz ▸ hxy, ih⟩

/-- `push_cigar` pushes every run of a CIGAR string, in order. -/
def pushAll (ops : Rle) (runs : Rle) : Rle := runs.foldl (fun acc r => pushOp r.1 r.2 acc) ops

theorem pushAll_spec : ∀ (runs ops : Rle), Canon ops →
    Canon (pushAll ops runs) ∧ expand (pushAll ops runs) = expand ops ++ expand runs
  | [], ops, h => ⟨h, by simp [pushAll, expand]⟩
  | r :: runs, ops, h => by
    obtain ⟨h1, h2⟩ := pushAll_spec runs (pushOp r.1 r.2 ops) (canon_pushOp _ _ _ h)
    refine ⟨h1, ?_⟩
    simp only [pushAll, List.foldl_cons] at h1 h2 ⊢
    rw [h2, expand_pushOp]
    simp [expand]

/-! ## The NM a CIGAR implies -/

variable {α : Type*} (mm : α → α → ℕ)

/-- `cigar_ed_implied` column by column, with `none` where the code returns None for an op other than M, I or D or a run past either string. -/
def edCols : List COp → List α → List α → Option ℕ
  | [], _, _ => some 0
  | .M :: ops, a :: q, b :: t => (mm a b + ·) <$> edCols ops q t
  | .I :: ops, _ :: q, t => (1 + ·) <$> edCols ops q t
  | .D :: ops, q, _ :: t => (1 + ·) <$> edCols ops q t
  | _, _, _ => none

/-- `cigar_ed_implied`, which also rejects a run of length 0. -/
def edImplied (c : Rle) (q t : List α) : Option ℕ :=
  if c.any (·.1 = 0) then none else edCols mm (expand c) q t

/-- The M/I/D columns, in `Align`'s type. -/
def toAlign : COp → Option Align.Op
  | .M => some .M
  | .I => some .I
  | .D => some .D
  | _ => none

theorem edCols_eq_cost : ∀ (ops : List COp) (q t : List α) (e : ℕ), edCols mm ops q t = some e →
    e = Align.cost mm (ops.filterMap toAlign) q t
  | [], _, _, e, h => by simp [edCols] at h; simp [Align.cost, h]
  | .M :: ops, a :: q, b :: t, e, h => by
    simp only [edCols, Option.map_eq_map, Option.map_eq_some_iff] at h
    obtain ⟨e', he', rfl⟩ := h
    simp [toAlign, Align.cost, edCols_eq_cost ops q t e' he']
  | .I :: ops, _ :: q, t, e, h => by
    simp only [edCols, Option.map_eq_map, Option.map_eq_some_iff] at h
    obtain ⟨e', he', rfl⟩ := h
    simp [toAlign, Align.cost, edCols_eq_cost ops q t e' he']
  | .D :: ops, [], _ :: t, e, h => by
    simp only [edCols, Option.map_eq_map, Option.map_eq_some_iff] at h
    obtain ⟨e', he', rfl⟩ := h
    simp [toAlign, Align.cost, edCols_eq_cost ops [] t e' he']
  | .D :: ops, a :: q, _ :: t, e, h => by
    simp only [edCols, Option.map_eq_map, Option.map_eq_some_iff] at h
    obtain ⟨e', he', rfl⟩ := h
    simp [toAlign, Align.cost, edCols_eq_cost ops (a :: q) t e' he']
  | .M :: _, [], _, _, h | .M :: _, _ :: _, [], _, h | .I :: _, [], _, _, h | .D :: _, _, [], _, h
  | .N :: _, _, _, _, h | .S :: _, _, _, _, h | .H :: _, _, _, _, h | .Eq :: _, _, _, _, h
  | .X :: _, _, _, _, h => by simp [edCols] at h

/-- The NM recomputed from a CIGAR is the edit cost of its columns, which the left-alignment proofs show never grows. -/
theorem edImplied_eq_cost {c : Rle} {q t : List α} {e : ℕ} (h : edImplied mm c q t = some e) :
    e = Align.cost mm ((expand c).filterMap toAlign) q t := by
  unfold edImplied at h
  split_ifs at h
  exact edCols_eq_cost mm _ q t e h

/-! ## clip_ends -/

/-- The columns `clip_ends` builds, each an op and whether it is a match, or `none` where it returns None. -/
def colsOf (isMatch : α → α → Bool) : List COp → List α → List α → Option (List (COp × Bool))
  | [], _, _ => some []
  | .M :: ops, a :: q, b :: t => ((.M, isMatch a b) :: ·) <$> colsOf isMatch ops q t
  | .I :: ops, _ :: q, t => ((.I, false) :: ·) <$> colsOf isMatch ops q t
  | .D :: ops, q, _ :: t => ((.D, false) :: ·) <$> colsOf isMatch ops q t
  | _, _, _ => none

/-- The query and text bases a list of columns consumes. -/
def qcols (cs : List (COp × Bool)) : ℕ := cs.countP (·.1 ≠ .D)
def tcols (cs : List (COp × Bool)) : ℕ := cs.countP (·.1 ≠ .I)

/-- The edit distance `clip_ends` counts, one per mismatch and one per gap base. -/
def nm (cs : List (COp × Bool)) : ℕ := cs.countP fun c => c.1 ≠ .M || !c.2

variable (isMatch : α → α → Bool)

/-- Built columns are the CIGAR's columns, and they consume exactly the bases they pair. -/
theorem colsOf_split : ∀ (ops : List COp) (q t : List α) (cs : List (COp × Bool)), colsOf isMatch ops q t = some cs →
    ∀ i, colsOf isMatch (ops.drop i) (q.drop (qcols (cs.take i))) (t.drop (tcols (cs.take i))) = some (cs.drop i)
  | ops, q, t, cs, h, 0 => by simpa [qcols, tcols] using h
  | [], _, _, cs, h, i + 1 => by simp [colsOf] at h; subst h; simp [colsOf]
  | .M :: ops, a :: q, b :: t, cs, h, i + 1 => by
    simp only [colsOf, Option.map_eq_map, Option.map_eq_some_iff] at h
    obtain ⟨cs', h', rfl⟩ := h
    simpa [qcols, tcols, List.countP_cons] using colsOf_split ops q t cs' h' i
  | .I :: ops, _ :: q, t, cs, h, i + 1 => by
    simp only [colsOf, Option.map_eq_map, Option.map_eq_some_iff] at h
    obtain ⟨cs', h', rfl⟩ := h
    simpa [qcols, tcols, List.countP_cons] using colsOf_split ops q t cs' h' i
  | .D :: ops, q, _ :: t, cs, h, i + 1 => by
    simp only [colsOf, Option.map_eq_map, Option.map_eq_some_iff] at h
    obtain ⟨cs', h', rfl⟩ := h
    have := colsOf_split ops q t cs' h' i
    simp only [List.take_succ_cons, List.drop_succ_cons, qcols, tcols, List.countP_cons] at this ⊢
    simpa using this
  | .M :: _, [], _, _, h, _ + 1 | .M :: _, _ :: _, [], _, h, _ + 1 | .I :: _, [], _, _, h, _ + 1
  | .D :: _, _, [], _, h, _ + 1 | .N :: _, _, _, _, h, _ + 1 | .S :: _, _, _, _, h, _ + 1
  | .H :: _, _, _, _, h, _ + 1 | .Eq :: _, _, _, _, h, _ + 1 | .X :: _, _, _, _, h, _ + 1 => by simp [colsOf] at h

/-- With the code's 0/1 mismatch cost, the columns' NM is the edit cost `edCols` gives. -/
theorem nm_eq_edCols : ∀ (ops : List COp) (q t : List α) (cs : List (COp × Bool)), colsOf isMatch ops q t = some cs →
    edCols (fun a b => if isMatch a b then 0 else 1) ops q t = some (nm cs)
  | [], _, _, cs, h => by simp [colsOf] at h; subst h; simp [edCols, nm]
  | .M :: ops, a :: q, b :: t, cs, h => by
    simp only [colsOf, Option.map_eq_map, Option.map_eq_some_iff] at h
    obtain ⟨cs', h', rfl⟩ := h
    simp only [edCols, nm_eq_edCols ops q t cs' h', Option.map_eq_map, Option.map_some, nm, List.countP_cons]
    cases isMatch a b <;> simp [Nat.add_comm]
  | .I :: ops, _ :: q, t, cs, h => by
    simp only [colsOf, Option.map_eq_map, Option.map_eq_some_iff] at h
    obtain ⟨cs', h', rfl⟩ := h
    simp only [edCols, nm_eq_edCols ops q t cs' h', Option.map_eq_map, Option.map_some, nm, List.countP_cons]
    simp; omega
  | .D :: ops, [], _ :: t, cs, h => by
    simp only [colsOf, Option.map_eq_map, Option.map_eq_some_iff] at h
    obtain ⟨cs', h', rfl⟩ := h
    simp only [edCols, nm_eq_edCols ops [] t cs' h', Option.map_eq_map, Option.map_some, nm, List.countP_cons]
    simp; omega
  | .D :: ops, a :: q, _ :: t, cs, h => by
    simp only [colsOf, Option.map_eq_map, Option.map_eq_some_iff] at h
    obtain ⟨cs', h', rfl⟩ := h
    simp only [edCols, nm_eq_edCols ops (a :: q) t cs' h', Option.map_eq_map, Option.map_some, nm, List.countP_cons]
    simp; omega
  | .M :: _, [], _, _, h | .M :: _, _ :: _, [], _, h | .I :: _, [], _, _, h | .D :: _, _, [], _, h
  | .N :: _, _, _, _, h | .S :: _, _, _, _, h | .H :: _, _, _, _, h | .Eq :: _, _, _, _, h
  | .X :: _, _, _, _, h => by simp [colsOf] at h

theorem colsOf_length : ∀ (ops : List COp) (q t : List α) (cs : List (COp × Bool)), colsOf isMatch ops q t = some cs →
    cs.length = ops.length ∧ qcols cs ≤ q.length
  | [], _, _, cs, h => by simp [colsOf] at h; subst h; simp [qcols]
  | .M :: ops, a :: q, b :: t, cs, h => by
    simp only [colsOf, Option.map_eq_map, Option.map_eq_some_iff] at h
    obtain ⟨cs', h', rfl⟩ := h
    have := colsOf_length ops q t cs' h'
    simp [qcols] at this ⊢; omega
  | .I :: ops, _ :: q, t, cs, h => by
    simp only [colsOf, Option.map_eq_map, Option.map_eq_some_iff] at h
    obtain ⟨cs', h', rfl⟩ := h
    have := colsOf_length ops q t cs' h'
    simp [qcols] at this ⊢; omega
  | .D :: ops, q, _ :: t, cs, h => by
    simp only [colsOf, Option.map_eq_map, Option.map_eq_some_iff] at h
    obtain ⟨cs', h', rfl⟩ := h
    have := colsOf_length ops q t cs' h'
    simp [qcols] at this ⊢; omega
  | .M :: _, [], _, _, h | .M :: _, _ :: _, [], _, h | .I :: _, [], _, _, h | .D :: _, _, [], _, h
  | .N :: _, _, _, _, h | .S :: _, _, _, _, h | .H :: _, _, _, _, h | .Eq :: _, _, _, _, h
  | .X :: _, _, _, _, h => by simp [colsOf] at h

/-- `clip_ends` with cut points `lcut` and `rcut`, whichever ones the scores chose.
It returns the S-clipped CIGAR's left clip, core columns and right clip, the new offset and the core's NM. -/
def clip (cs : List (COp × Bool)) (toff lcut rcut : ℕ) : ℕ × List (COp × Bool) × ℕ × ℕ × ℕ :=
  let core := (cs.drop lcut).take (cs.length - lcut - rcut)
  (qcols (cs.take lcut), core, qcols (cs.drop (cs.length - rcut)), toff + tcols (cs.take lcut), nm core)

theorem qcols_append (a b : List (COp × Bool)) : qcols (a ++ b) = qcols a + qcols b := by simp [qcols]

/-- The left clip, core and right clip together cover the whole read, so SEQ and the CIGAR still agree. -/
theorem clip_query_length (cs : List (COp × Bool)) (toff lcut rcut : ℕ) (h : lcut + rcut < cs.length) :
    let r := clip cs toff lcut rcut
    r.1 + qcols r.2.1 + r.2.2.1 = qcols cs := by
  simp only [clip]
  have e : cs = cs.take lcut ++ ((cs.drop lcut).take (cs.length - lcut - rcut) ++ cs.drop (cs.length - rcut)) := by
    conv_lhs => rw [← List.take_append_drop lcut cs]
    congr 1
    conv_lhs => rw [← List.take_append_drop (cs.length - lcut - rcut) (cs.drop lcut)]
    congr 1
    rw [List.drop_drop]; congr 1; omega
  conv_rhs => rw [e]
  rw [qcols_append, qcols_append, Nat.add_assoc]

/-- The columns of a prefix of the ops are the prefix of the columns. -/
theorem colsOf_take : ∀ (ops : List COp) (q t : List α) (cs : List (COp × Bool)), colsOf isMatch ops q t = some cs →
    ∀ i, colsOf isMatch (ops.take i) q t = some (cs.take i)
  | _, _, _, _, _, 0 => by simp [colsOf]
  | [], _, _, cs, h, i + 1 => by simp [colsOf] at h; subst h; simp [colsOf]
  | .M :: ops, a :: q, b :: t, cs, h, i + 1 => by
    simp only [colsOf, Option.map_eq_map, Option.map_eq_some_iff] at h
    obtain ⟨cs', h', rfl⟩ := h
    simp [colsOf, colsOf_take ops q t cs' h' i]
  | .I :: ops, _ :: q, t, cs, h, i + 1 => by
    simp only [colsOf, Option.map_eq_map, Option.map_eq_some_iff] at h
    obtain ⟨cs', h', rfl⟩ := h
    simp [colsOf, colsOf_take ops q t cs' h' i]
  | .D :: ops, [], _ :: t, cs, h, i + 1 => by
    simp only [colsOf, Option.map_eq_map, Option.map_eq_some_iff] at h
    obtain ⟨cs', h', rfl⟩ := h
    simp [colsOf, colsOf_take ops [] t cs' h' i]
  | .D :: ops, a :: q, _ :: t, cs, h, i + 1 => by
    simp only [colsOf, Option.map_eq_map, Option.map_eq_some_iff] at h
    obtain ⟨cs', h', rfl⟩ := h
    simp [colsOf, colsOf_take ops (a :: q) t cs' h' i]
  | .M :: _, [], _, _, h, _ + 1 | .M :: _, _ :: _, [], _, h, _ + 1 | .I :: _, [], _, _, h, _ + 1
  | .D :: _, _, [], _, h, _ + 1 | .N :: _, _, _, _, h, _ + 1 | .S :: _, _, _, _, h, _ + 1
  | .H :: _, _, _, _, h, _ + 1 | .Eq :: _, _, _, _, h, _ + 1 | .X :: _, _, _, _, h, _ + 1 => by simp [colsOf] at h

/-- The NM `clip_ends` reports is the edit cost of the kept core's ops, against the query after the left clip and the text from the new offset. -/
theorem clip_nm {ops : List COp} {q t : List α} {cs : List (COp × Bool)} {toff : ℕ} (lcut rcut : ℕ)
    (h : colsOf isMatch ops q (t.drop toff) = some cs) :
    let r := clip cs toff lcut rcut
    edCols (fun a b => if isMatch a b then 0 else 1) ((ops.drop lcut).take (cs.length - lcut - rcut))
      (q.drop r.1) (t.drop r.2.2.2.1) = some r.2.2.2.2 := by
  simp only [clip]
  have hs := colsOf_split isMatch ops q _ cs h lcut
  rw [List.drop_drop] at hs
  exact nm_eq_edCols isMatch _ _ _ _ (colsOf_take isMatch _ _ _ _ hs _)

/-- Clipping never raises NM. -/
theorem clip_nm_le (cs : List (COp × Bool)) (toff lcut rcut : ℕ) : (clip cs toff lcut rcut).2.2.2.2 ≤ nm cs := by
  simp only [clip, nm]
  exact List.Sublist.countP_le ((List.take_sublist _ _).trans (List.drop_sublist _ _))

/-! ## hard_clip -/

/-- `hard_clip` turns soft clips of `lead` and `trail` bases around `body` into hard clips. -/
def hardClip (lead trail : ℕ) (body : Rle) : Rle := (lead, .H) :: body ++ [(trail, .H)]

/-- The hard-clipped CIGAR consumes exactly the bases left in the cut sequence. -/
theorem hardClip_qlen {lead trail : ℕ} {body : Rle} {seq : List α}
    (h : qlen ((lead, .S) :: body ++ [(trail, .S)]) = seq.length) :
    qlen (hardClip lead trail body) = ((seq.drop lead).take (seq.length - lead - trail)).length := by
  simp only [qlen, hardClip, List.cons_append, runSum, runSum_append, COp.query] at h ⊢
  simp only [List.length_take, List.length_drop]
  simp at h ⊢
  omega

/-- Soft and hard clips both consume no reference, so POS, END and TLEN do not change. -/
theorem hardClip_refSpan (lead trail : ℕ) (body : Rle) :
    refSpan (hardClip lead trail body) = refSpan ((lead, .S) :: body ++ [(trail, .S)]) := by
  simp [refSpan, hardClip, runSum, runSum_append, COp.ref]

end OptikProofs.Cigar
