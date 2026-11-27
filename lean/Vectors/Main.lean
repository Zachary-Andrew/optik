import OptikProofs.Spss
import OptikProofs.Dbg
import OptikProofs.Cigar
import Lean.Data.Json

/-!
# Test vectors from the Lean models

`lake exe optik-vectors --out DIR [--seed N]` writes `spss.json`, `dbg.json`, `dbg2.json`, `mapq.json`, `cigar.json` and `bpm.json` into DIR.
`bench/validate_lean.py` runs the optik binary on the first three, and Rust unit tests check the others against `posterior_mapq_w`, the CIGAR helpers and the bit-parallel placement.

Expected answers come from the same definitions the theorems are about.
Spss presence is `Spss.occursB`, proved equal to `Spss.Occurs`, so by `Spss.present_iff` it is exactly what lookup must answer.
Dbg unitigs come from a walker that decides `Dbg.Compacts` (Definition 7.1) over the edges the reads contain.
The walker itself is executable code rather than a proved function, so the Dbg vectors are a second implementation of the definition and not a proof about it.
At -c 2 (`dbg2.json`) the walker uses only edges between solid k-mers, `Count.TrimEdge`, and the vectors also list those edges so the GFA links can be checked.
MAPQ is computed in floating point, step for step as `posterior_mapq_w` computes it.
CIGAR answers are `Cigar.edImplied`, `Cigar.refSpan` and `Cigar.colCount`, so the Rust test also checks the code's digit parsing.
Placement answers come from the recurrence of `Bpm.D`, filled one column at a time, which is what `Bpm.bpm_eq_dp` says the code returns.
-/

open Lean OptikProofs OptikProofs.Kmer

/-- A fixed-seed linear congruential generator, so the same seed always gives the same vectors. -/
structure Rng where
  state : UInt64

def Rng.next (r : Rng) : Rng × Nat :=
  let s := r.state * 6364136223846793005 + 1442695040888963407
  (⟨s⟩, (s >>> 33).toNat)

def Rng.below (r : Rng) (n : Nat) : Rng × Nat :=
  let (r, v) := r.next
  (r, if n = 0 then 0 else v % n)

def Rng.bases (r : Rng) (len : Nat) : Rng × List Base := Id.run do
  let mut r := r
  let mut out := []
  for _ in [0:len] do
    let (r', v) := r.below 4
    r := r'
    out := ⟨v % 4, Nat.mod_lt _ (by decide)⟩ :: out
  return (r, out)

def baseChar (b : Base) : Char := "ACGT".toList.getD b.val 'N'

def str (x : List Base) : String := String.ofList (x.map baseChar)

/-! ## Spss -/

def spssCase (r : Rng) : Rng × Json := Id.run do
  let mut r := r
  let (r1, k) := r.below 20
  r := r1
  let k := k + 11
  let (r2, nstr) := r.below 4
  r := r2
  let mut strs : List (List Base) := []
  for _ in [0:nstr + 1] do
    let (r3, len) := r.below 60
    let (r4, s) := r3.bases (len + k)
    r := r4
    strs := strs ++ [s]
  let mut queries : List Json := []
  for i in [0:12] do
    let (r5, which) := r.below 3
    let (r6, sid) := r5.below strs.length
    let s := strs.getD sid []
    let (r7, off) := r6.below (s.length + 1 - k)
    let (r8, rnd) := r7.bases k
    r := r8
    let w := Spss.window k s off
    let x := if which = 0 then w else if which = 1 then rc w else rnd
    let _ := i
    let expect := Spss.occursB k strs x || Spss.occursB k strs (rc x)
    queries := queries ++ [Json.mkObj [("kmer", str x), ("present", toJson expect)]]
  return (r, Json.mkObj [("k", toJson k), ("strings", toJson (strs.map str)), ("queries", Json.arr queries.toArray)])

/-! ## Dbg -/

def allBases : List Base := [0, 1, 2, 3]

/-- The edge from x to X exists when the reads contain its (k+1)-mer on either strand, as optik's symmetric masks record it. -/
def hasEdge (k : Nat) (words : List (List Base)) (x X : List Base) : Bool :=
  let w := Dbg.rightWord k x X
  x.tail == X.dropLast && (words.contains w || words.contains (rc w))

/-- A k-mer is solid when its canonical form occurs at least `minCount` times in the reads, on either strand. -/
def isSolid (minCount : Nat) (kmers : List (List Base)) (x : List Base) : Bool :=
  (kmers.filter fun y => canon y == canon x).length ≥ minCount

/-- An observed edge whose two k-mers are both solid, as `Count.TrimEdge` defines it. -/
def solidEdge (k : Nat) (words : List (List Base)) (solid : List Base → Bool) (x X : List Base) : Bool :=
  hasEdge k words x X && solid x && solid X

def succs (k : Nat) (words : List (List Base)) (solid : List Base → Bool) (x : List Base) : List (List Base) :=
  (allBases.map fun b => x.tail ++ [b]).filter (solidEdge k words solid x)

def preds (k : Nat) (words : List (List Base)) (solid : List Base → Bool) (X : List Base) : List (List Base) :=
  (allBases.map fun b => b :: X.dropLast).filter (fun z => solidEdge k words solid z X)

/-- Definition 7.1, decided over the edges in `words` between solid k-mers. -/
def compacts (k : Nat) (words : List (List Base)) (solid : List Base → Bool) (x X : List Base) : Bool :=
  succs k words solid x == [X] && preds k words solid X == [x] && canon X != canon x && x != rc x && X != rc X

def canonForm (x : List Base) : List Base := if pack x ≤ pack (rc x) then x else rc x

/-- Follows compaction steps from z, at most `fuel` of them, returning the canonical k-mers visited. -/
def walk (k : Nat) (words : List (List Base)) (solid : List Base → Bool) (z : List Base) : Nat → List (List Base)
  | 0 => [canonForm z]
  | fuel + 1 =>
    match succs k words solid z with
    | [X] => if compacts k words solid z X && canonForm X != canonForm z then canonForm z :: walk k words solid X fuel
             else [canonForm z]
    | _ => [canonForm z]

def windowsOf (k : Nat) (s : List Base) : List (List Base) :=
  (List.range (s.length + 1 - k)).map fun i => Spss.window k s i

/-- The unitigs over the k-mers that occur at least `minCount` times. -/
def dbgUnitigs (k : Nat) (reads : List (List Base)) (minCount : Nat := 1) : List (List String) := Id.run do
  let kmers := reads.flatMap (windowsOf k)
  let solid := isSolid minCount kmers
  let nodes := (kmers.map canonForm).eraseDups.filter solid
  let words := reads.flatMap (windowsOf (k + 1))
  let n := nodes.length
  let mut covered : List (List Base) := []
  let mut unitigs : List (List String) := []
  -- Starts are oriented k-mers nothing compacts into, and every non-cycle unitig is walked from one of them.
  for c in nodes do
    for z in [c, rc c] do
      if (preds k words solid z).all (fun w => !compacts k words solid w z) && !covered.contains (canonForm z) then
        let path := (walk k words solid z n).eraseDups
        covered := covered ++ path
        unitigs := unitigs ++ [(path.map str).mergeSort (· ≤ ·)]
  -- What remains lies on pure cycles, each walked from any of its k-mers.
  for c in nodes do
    if !covered.contains c then
      let path := (walk k words solid c n).eraseDups.filter (fun y => !covered.contains y)
      covered := covered ++ path
      unitigs := unitigs ++ [(path.map str).mergeSort (· ≤ ·)]
  return unitigs.mergeSort (fun a b => a.headD "" ≤ b.headD "")

/-- A random k and read set for the Dbg vectors. -/
def dbgReads (r : Rng) : Rng × Nat × List (List Base) := Id.run do
  let mut r := r
  let (r1, kk) := r.below 6
  r := r1
  -- k from 4 to 9, so even k (and with it palindromic k-mers) comes up often.
  let k := kk + 4
  let (r2, seedLen) := r.below 30
  let (r3, genome) := r2.bases (seedLen + 3 * k)
  r := r3
  let mut reads : List (List Base) := [genome]
  -- Overlapping pieces and a reverse-complemented copy add branches, repeats and hairpins.
  for _ in [0:3] do
    let (r4, off) := r.below (genome.length + 1 - k)
    let (r5, len) := r4.below 12
    let (r6, flip) := r5.below 2
    r := r6
    let piece := (genome.drop off).take (len + k)
    reads := reads ++ [if flip = 0 then piece else rc piece]
  let (r7, extra) := r.bases (k + 6)
  r := r7
  reads := reads ++ [extra ++ rc extra]
  return (r, k, reads)

def dbgCase (r : Rng) : Rng × Json :=
  let (r, k, reads) := dbgReads r
  (r, Json.mkObj [("k", toJson k), ("reads", toJson (reads.map str)), ("unitigs", toJson (dbgUnitigs k reads))])

/-- The canonical (k+1)-mers of the edges between solid k-mers, each as the smaller of the word and its reverse complement. -/
def solidEdgeWords (k : Nat) (reads : List (List Base)) (minCount : Nat) : List String :=
  let solid := isSolid minCount (reads.flatMap (windowsOf k))
  let words := (reads.flatMap (windowsOf (k + 1))).filter fun w => solid (w.take k) && solid (w.drop 1)
  ((words.map fun w => str (canonForm w)).eraseDups).mergeSort (· ≤ ·)

/-- A `dbgCase` read set with two copies of part of it and one of them with a substitution, counted at -c 2.
The doubled part gives solid k-mers, and the substitution gives k-mers that occur once, whose edges the trim must drop. -/
def dbg2Case (r : Rng) : Rng × Json := Id.run do
  let (r, k, reads0) := dbgReads r
  let genome := reads0.headD []
  let mut r := r
  let (r1, off) := r.below (genome.length + 1 - k)
  let (r2, pos) := r1.below (2 * k)
  let (r3, sub) := r2.below 3
  r := r3
  let piece := (genome.drop off).take (3 * k)
  let old := piece.getD pos 0
  let mutated := piece.set pos ⟨(old.val + sub + 1) % 4, Nat.mod_lt _ (by decide)⟩
  let reads := reads0 ++ [piece, mutated]
  return (r, Json.mkObj [("k", toJson k), ("min_count", toJson 2), ("reads", toJson (reads.map str)),
    ("unitigs", toJson (dbgUnitigs k reads 2)), ("edges", toJson (solidEdgeWords k reads 2))])

/-! ## MAPQ -/

/-- `edit_ratio` in tools/src/map/paired.rs, with its fitted constants. -/
def editRatio (edits : Int) (len : Nat) : Float :=
  let e := (Float.ofInt edits + 1.2608) / Float.ofNat (max len 1) + 0.016
  let r := 0.8228 * e / max (1.0 - e) 1e-3
  min (max r 1e-6) 0.9

def powi (x : Float) (n : Int) : Float :=
  if n ≥ 0 then x ^ Float.ofInt n else 1.0 / x ^ Float.ofInt (-n)

/-- `posterior_mapq_w` with its unrounded value, so cases next to a rounding boundary can be skipped. -/
def mapqRaw (best : Int) (comps : List (Int × Float)) (len : Nat) (hidden : Float) : Option Float :=
  if comps.any (fun c => c.1 ≤ best && c.2 ≥ 1.0) then none else
  let ratio := editRatio best len
  let odds := (comps.map fun c => c.2 * powi ratio (c.1 - best)).foldl (· + ·) 0.0 + hidden
  let p := max (odds / (1.0 + odds)) 1e-9
  some (-10.0 * Float.log10 p)

def mapqOf (raw : Option Float) : Int :=
  match raw with
  | none => 0
  | some v => max 0 (min 60 (Float.round v).toInt64.toInt)

def nearHalf (raw : Option Float) : Bool :=
  match raw with
  | none => false
  | some v => (v - v.floor - 0.5).abs < 1e-6

def mapqCase (r : Rng) : Rng × Option Json := Id.run do
  let mut r := r
  let (r1, best) := r.below 12
  let (r2, ncomp) := r1.below 4
  let (r3, len) := r2.below 300
  let (r4, hid) := r3.below 1000
  r := r4
  let mut comps : List (Int × Float) := []
  for _ in [0:ncomp] do
    let (r5, d) := r.below 14
    let (r6, wq) := r5.below 5
    r := r6
    comps := comps ++ [((d : Int), Float.ofNat wq * 0.25 + 0.25)]
  let len := len + 50
  let hidden := Float.ofNat hid * 1e-4
  let raw := mapqRaw best comps len hidden
  if nearHalf raw then return (r, none)
  return (r, some (Json.mkObj [("best_ed", toJson (best : Int)),
    ("competitors", Json.arr (comps.map fun c => Json.arr #[toJson c.1, toJson c.2]).toArray),
    ("len", toJson len), ("hidden", toJson hidden), ("mapq", toJson (mapqOf raw))]))

/-! ## CIGAR -/

open Cigar in
def opChar : COp → Char
  | .M => 'M' | .I => 'I' | .D => 'D' | .N => 'N' | .S => 'S' | .H => 'H' | .Eq => '=' | .X => 'X'

/-- Under the code's mismatch rule a base pair costs 0 only when both bases agree and are ACGT. -/
def ntCost (a b : Char) : Nat := if a == b && "ACGT".toList.contains b then 0 else 1

def rleString (c : Cigar.Rle) : String := String.join (c.map fun r => s!"{r.1}{opChar r.2}")

/-- A random M/I/D CIGAR with a query and text it mostly fits.
Some cases have a zero-length run, an S op or a text one base short, where the code must answer None. -/
def cigarCase (r : Rng) : Rng × Json := Id.run do
  let mut r := r
  let (r1, nruns) := r.below 6
  r := r1
  let mut runs : Cigar.Rle := []
  for _ in [0:nruns + 1] do
    let (r2, op) := r.below 3
    let (r3, len) := r2.below 12
    r := r3
    runs := runs ++ [(len + 1, if op = 0 then .M else if op = 1 then .I else .D)]
  let (r4, fault) := r.below 10
  let (r5, toff) := r4.below 4
  r := r5
  if fault = 0 then runs := runs ++ [(0, .M)]
  if fault = 1 then runs := runs ++ [(3, .S)]
  -- The query and text follow the columns, with a text base copying the query base three times in four.
  let mut q : List Char := []
  let mut t : List Char := []
  for (n, o) in runs do
    for _ in [0:n] do
      let (r6, b) := r.below 5
      let (r7, keep) := r6.below 4
      let (r8, b') := r7.below 4
      r := r8
      let qb := "ACGTN".toList.getD b 'N'
      let tb := if keep = 0 then "ACGT".toList.getD b' 'A' else qb
      if o = .M then q := q ++ [qb]; t := t ++ [tb]
      else if o = .I then q := q ++ [qb]
      else if o = .D then t := t ++ [tb]
      else if o = .S then q := q ++ [qb]
  let (r9, pad) := r.bases toff
  r := r9
  let text := (pad.map baseChar) ++ (if fault = 2 then t.dropLast else t)
  let nm := Cigar.edImplied ntCost runs q (text.drop toff)
  return (r, Json.mkObj [("cigar", toJson (rleString runs)), ("query", toJson (String.ofList q)),
    ("text", toJson (String.ofList text)), ("toff", toJson toff), ("nm", toJson nm),
    ("ref_span", toJson (Cigar.refSpan runs)), ("cols", toJson (Cigar.colCount runs))])

/-! ## Bit-parallel placement -/

/-- The code's base codes: A, C, G and T in either case match only themselves, and anything else matches nothing. -/
def bpmCost (a b : Char) : Int :=
  if a.toUpper == b.toUpper && "ACGT".toList.contains b.toUpper then 0 else 1

/-- `Bpm.D j m` for every column j, by the DP's recurrence one column at a time. -/
def bottomRow (q t : List Char) : List Int := Id.run do
  let m := q.length
  let mut col : Array Int := (Array.range (m + 1)).map Int.ofNat
  let mut out : List Int := [col[m]!]
  for b in t do
    let mut next : Array Int := #[0]
    for r in [0:m] do
      let c := bpmCost (q.getD r 'N') b
      next := next.push (min (min (col[r]! + c) (col[r + 1]! + 1)) (next[r]! + 1))
    col := next
    out := out ++ [col[m]!]
  return out

def bpmChar (r : Rng) : Rng × Char :=
  let (r, v) := r.below 24
  (r, if v < 20 then "ACGT".toList.getD (v % 4) 'A' else if v < 22 then "acgt".toList.getD (v % 4) 'a' else 'N')

def bpmString (r : Rng) (len : Nat) : Rng × List Char := Id.run do
  let mut r := r
  let mut out := []
  for _ in [0:len] do
    let (r', c) := bpmChar r
    r := r'
    out := out ++ [c]
  return (r, out)

/-- A query of 1 to 200 bases, so one to four words, against one to four texts that mostly hold an edited copy of it. -/
def bpmCase (r : Rng) : Rng × Json := Id.run do
  let mut r := r
  let (r1, qlen) := r.below 200
  let (r2, ntexts) := r1.below 4
  let (r3, maxEd) := r2.below 40
  r := r3
  let (r4, q) := bpmString r (qlen + 1)
  r := r4
  let mut texts : List Json := []
  for _ in [0:ntexts + 1] do
    let (r5, left) := r.below 40
    let (r6, right) := r5.below 40
    let (r7, nedit) := r6.below 12
    let (r8, a) := bpmString r7 left
    let (r9, z) := bpmString r8 right
    r := r9
    let mut body := q
    for _ in [0:nedit] do
      let (r10, pos) := r.below (body.length + 1)
      let (r11, kind) := r10.below 3
      let (r12, c) := bpmChar r11
      r := r12
      body := if kind = 0 then body.set pos c else if kind = 1 then body.take pos ++ [c] ++ body.drop pos
        else body.take pos ++ body.drop (pos + 1)
    let t := a ++ body ++ z
    let row := bottomRow q t
    let best := row.foldl min (row.headD 0)
    -- The leftmost exclusive end with the least distance, as bpm_placements reports it.
    let ends := (row.drop 1).zipIdx.filter (fun p => p.1 == (row.drop 1).foldl min (row.getD 1 0))
    texts := texts ++ [Json.mkObj [("text", toJson (String.ofList t)), ("best", toJson (min best (maxEd + 1 : Int))),
      ("end", toJson ((ends.headD (0, 0)).2 + 1))]]
  return (r, Json.mkObj [("query", toJson (String.ofList q)), ("max_ed", toJson maxEd), ("texts", Json.arr texts.toArray)])

/-! ## Driver -/

def collect {α : Type} (r : Rng) (n : Nat) (gen : Rng → Rng × α) : Rng × List α := Id.run do
  let mut r := r
  let mut out := []
  for _ in [0:n] do
    let (r', v) := gen r
    r := r'
    out := out ++ [v]
  return (r, out)

def main (args : List String) : IO UInt32 := do
  let rec parse : List String → Option String → Nat → Option String × Nat
    | "--out" :: d :: rest, _, s => parse rest (some d) s
    | "--seed" :: v :: rest, o, _ => parse rest o (v.toNat?.getD 1)
    | _ :: rest, o, s => parse rest o s
    | [], o, s => (o, s)
  let (out?, seed) := parse args none 1
  let some out := out? | IO.eprintln "usage: optik-vectors --out DIR [--seed N]"; return 1
  IO.FS.createDirAll out
  let r : Rng := ⟨seed.toUInt64⟩
  let (r, spss) := collect r 40 spssCase
  let (r, dbg) := collect r 40 dbgCase
  let (r, mapq) := collect r 400 mapqCase
  let (r, cigar) := collect r 300 cigarCase
  let (r, dbg2) := collect r 40 dbg2Case
  let (_, bpm) := collect r 150 bpmCase
  IO.FS.writeFile (out ++ "/spss.json") (Json.arr spss.toArray).pretty
  IO.FS.writeFile (out ++ "/dbg.json") (Json.arr dbg.toArray).pretty
  IO.FS.writeFile (out ++ "/mapq.json") (Json.arr (mapq.filterMap id).toArray).pretty
  IO.FS.writeFile (out ++ "/cigar.json") (Json.arr cigar.toArray).pretty
  IO.FS.writeFile (out ++ "/dbg2.json") (Json.arr dbg2.toArray).pretty
  IO.FS.writeFile (out ++ "/bpm.json") (Json.arr bpm.toArray).pretty
  IO.println s!"wrote {spss.length} spss, {dbg.length} dbg, {(mapq.filterMap id).length} mapq, {cigar.length} cigar, {dbg2.length} dbg -c 2 and {bpm.length} bpm cases to {out}"
  return 0
