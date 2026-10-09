# Opt-in size-weighted adult production candidate (2026-10-08)

## Hypothesis and state contract

An outbreak can leave many counted adults but relatively few large breeders.
If newly matured adults (15–25 cm) produce fewer gametes than >25 cm adults,
the reproductive trough may deepen while a subsequent cohort grows and can
support a second peak. This is a testable candidate, **not** an assertion that
1991 or inter-wave size distributions are known. The 3 COTS/ha Allee
half-saturation remains fixed; tow conversion and observed peak detection do
not change. The model's `N[3]` stays total adults in COTS/ha, consistent with
the >=15 cm observation contract.

Track `L`, large adults (>25 cm) in COTS/ha, with small adults `S=N[3]-L`
(15–25 cm). For the pre-dispersal potential-larvae calculation replace only
the breeder count:

```
B = w_S*S + L
w_S = exp(beta*(D_S-D_L))
P = fecundity_scale * body_condition^2 * B
    * exp(-b_ricker*N[3]) * N[3]^2/(theta^2+N[3]^2)
```

`P` has the current pre-pelagic production units, COTS ha^-1 yr^-1.
`beta=0.0115 mm^-1` is the female gonad-mass diameter slope in Babcock,
Milton & Pratchett (2016, DOI:10.1007/s00227-016-3009-5). `D_S=200 mm` and
`D_L=300 mm` are **representative class diameters**, not observed class means;
they give `w_S=0.317`. The large class is normalized to 1 so this test does
not quietly raise the existing per-adult fecundity scale. This ratio is a
relative gonad-mass proxy, not a measured conversion to settled recruits.
The paper sampled substantially larger reproductive animals and its scatter
is material (female fit R²=0.55), so diameter representatives and `beta`
must be sensitivity inputs, not fixed biological truths. Retain the existing
Allee function of **total** adults to isolate size weighting from a new
fertilisation law.

Each year, after baseline adult survival `s=(1-m3)*food_survival`, let
`g(F+S_coral)` be the probability that a surviving small adult moves to the
large class. The proposed minimal food-dependent form is
`g=g_max*clamp((F+S_coral)/growth_cover,0,1)`. Then
`L_next=s*L + g*s*(N[3]-L)` and
`N[3]_next=s*N[3]+matured_subadults`; newly matured subadults join the small
class. The same survival and grazing per adult as the control are retained:
any change in coral comes through recruitment and abundance, not an extra
consumption or mortality effect. The control is the existing unweighted
model path with no size-state feedback. The switch defaults off and is the
rollback path.

## Bounded qualitative test

Use the archived 1991 Owen domain, source/sink orientation, sampled annual
hydrodynamic years, evidence-derived external outbreak boundary, separated
recruitment and 3 COTS/ha Allee setting. First reproduce an archived control
with identical initial state, coral, `q`, and forcing seed. The smallest
factorial is the control versus (A) size weighting at fixed `g=0.5/year` and
(B) size weighting with food-mediated growth (`g_max=0.5/year`,
`growth_cover=0.4` coral fraction). This distinguishes size weighting from
coral-conditioned replenishment of large breeders. Run a *bounded* sensitivity
over initial large-adult fractions 0.25 and 0.75, and class representatives
200/300 versus 200/350 mm. These values are diagnostic; historical size and
growth rates are unobserved. Freeze the same three-year-smoothed LTMP/EOTR
COTS/tow scoring and one coherent provisional `q` for all arms.

Log small and large adult density, effective breeder density, growth flux,
pre-pelagic production, local retention, internal/external settlement,
total adult density, modeled tow CPUE, and coral cover, all by reef and year.
Use the existing peak timing, height, prominence and inter-peak trough audit.
The candidate fails the qualitative gate if it merely reduces the first peak,
does not get a distinct second peak closer to observations, or creates deeper
troughs by degrading coral dynamics or nonnegative/mass-accounting checks.
No large optimizer or promotion follows a one-seed pass: repeat seeds,
conversion sensitivities, held-out reefs/years, and regional validation are
subsequent gates. Preserve failed evaluations and provenance hashes.

## Repository boundary

The ecological transition and pure unit tests belong in sibling `COTSMod.jl`;
the ADRIA adapter, output logging and bounded screen belong here. Do not
emulate a core change by silently changing the ADRIA adapter or observation
model. Review the uncommitted edits in both repositories together before
merging or branching; the ADRIA sandbox resolves its local `COTSMod.jl` path.

## First paired result and decision

The first bounded screen is archived at
`runs/20261008Tsize_weighted_pair020_v2/`. Its unweighted control exactly
replayed the previous pair-20 adult, coral and hydrodynamic-year trajectories
to the script's `1e-12` tolerance. All five arms completed and passed
small+large=total adult accounting; the independent `size_audit.json` records
input hashes and per-reef metrics. No arm meets the qualitative gate. The
worst two-wave reef trough ratio is 0.360 in the control and 0.359 in the
food-mediated size arm, both above the <=0.10 requirement. First peaks remain
several years late, the four-reef score still misses one of seven observed
peaks, and modeled coral remains too depleted relative to 9 m manta data
(subject to spatial-support differences). Weighting clearly reduces
effective breeders and local/internal recruitment, while external boundary
settlement stays fixed. The supported Allee parameter remains at 3 COTS/ha.

Decision: retain `size_weighted_fecundity=false` as the default and do not
promote this mechanism. There is no ecological rollback migration because
the switch is opt-in; disable it to recover the legacy transition. The next
discriminating factorial should test whether a bounded post-peak adult
survival/grazing response, separately from external boundary suppression,
can convert lower recruitment into a true adult crash without erasing the
second peak or worsening coral. That new mechanism requires its own protocol,
pure tests and paired control before optimization.
