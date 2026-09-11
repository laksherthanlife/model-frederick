"""Make oxidative stress dosable, so the GEM can price it instead of manufacturing it.

Yeast9 v9.0.2 has hydrogen peroxide in four compartments and superoxide in two, and no way
to put either there from outside. Superoxide has NO producing reaction at all, so both
superoxide dismutases are dead; peroxide has no extracellular species and no exchange, so
the only way to raise cytosolic H2O2 was to run polyamine oxidase -- which is why forcing
catalase ``r_0255`` cost growth for MANUFACTURING peroxide rather than for clearing a dose.

Installed on a copy of the model:

* ``h2o2_e`` with an exchange and a passive plasma-membrane diffusion, so a medium dose
  reaches the cytosol;
* a mitochondrial Qo-site electron leak feeding matrix superoxide, carrying complex III's
  own GPR, and the Yno1p NADPH oxidase feeding cytosolic superoxide, which between them
  give both SODs a source;
* TSA1 (YML028W), the dominant cytosolic peroxiredoxin, onto ``r_0550`` beside TSA2.

Every reaction is mass- and charge-balanced and every coefficient is electron conservation
or cited in :data:`CITATIONS`. Nothing is forced: each addition defaults to zero flux, so
installing the module leaves growth unchanged and the ODE layer sets the flux to be priced.

ALL THREE ROS ROUTES DEFAULT TO (0, 0), shut in both directions, and only a caller opens
one. A route the LP may use at its own discretion is not "no ROS": it is a free sink or a
free bypass, and each of the three was measured doing exactly that.

The exchange at (0, 1000) was a free unlimited peroxide sink that rescued the
clearance-null mutant to wild type. The two superoxide sources at (0, 1000) were an
alternative terminal oxidase -- leak, then SOD, then cytochrome c peroxidase r_0437 --
which S. cerevisiae does not have. Borrowing complex III's GPR does NOT close that: a GPR
constrains no flux in a model without enzyme costs, and the bypass runs around complex IV,
not complex III. Measured at glucose 10, oxygen free: a COX1 or COX2 null grew 0.19130 ->
0.37194 /h and doubled its maximum ATP turnover 20 -> 40, and with every exchange open the
WILD TYPE gained 4.16% growth (77.556 -> 80.780 /h) on a leak of 656 mmol/gDCW/h. The GPR
is kept because it is true, but the zero upper bound is what makes the module safe.

THE PRICE LIST, all at glucose <= 10 with oxygen free (base growth 0.887685 /h):

    imposed flux                        growth cost           linear to     dies at
    H2O2 uptake, catalase present       0.0000 /h per mmol    500 mmol      never
    H2O2 uptake, catalase deleted       0.0074244 /h per mmol 50 mmol       114.053
    mitochondrial leak ROSLEAK_Q_m      0.0117542 /h per mmol 40 mmol       60.0
    NADPH oxidase NOX_YNO1_c            0.0074244 /h per mmol 50 mmol       114.053

The two thiol routes and the NADPH oxidase share a slope because they share a currency: one
NADPH per H2O2. 114.053 is the network's whole redox throughput at glucose 10, which is why
the peroxide ceiling and the oxidase ceiling are the same number.

WHAT THIS DOES NOT BUY. Catalase is cofactor-free, so with it present the GEM prices a dose
at exactly nothing. The measured cost appears only through the thiol routes. Peroxide
toxicity is Fenton chemistry and macromolecular damage, and a stoichiometric model has none
of it.

AND ONE TRAP, MEASURED. Catalase turns 2 H2O2 into 2 H2O and an O2, so in an OXYGEN-LIMITED
regime a dose is an oxygen supply that walks past the oxygen bound: at glucose 11.1 /
oxygen 3.7 a dose of 5 RAISES growth 22.31% and a dose of 50 raises it 185%, at which point
the cell SECRETES oxygen. Every cost above was measured with oxygen free for that reason.
The ODE layer must not read a dosed growth rate as toxicity where oxygen is binding.
"""

from __future__ import annotations

from dataclasses import dataclass

import cobra
import pandas as pd

from .solver import FVA_PROCESSES, growth_or_none

__all__ = [
    "CITATIONS",
    "CLEARANCE_REACTIONS",
    "CLEARANCE_ROUTES",
    "COMPLEX_III",
    "CYTOSOLIC_CATALASE",
    "MITOCHONDRIAL_LEAK",
    "NADPH_OXIDASE",
    "PEROXIDE_EXCHANGE",
    "PEROXIDE_TRANSPORT",
    "REDOX_THROUGHPUT_CEILING",
    "ROS_REACTION_IDS",
    "THIOREDOXIN_PEROXIDASE",
    "TSA1_ORF",
    "ClearanceCost",
    "add_ros_module",
    "clearance_cost",
    "clearance_ladder",
    "closed_medium_ranges",
    "dose_peroxide",
    "dose_superoxide",
    "energy_from_peroxide",
    "install_tsa1",
    "peroxide_dose_response",
    "superoxide_leak_response",
]

PEROXIDE_EXCHANGE = "EX_h2o2_e"
PEROXIDE_TRANSPORT = "H2O2t_pm"
MITOCHONDRIAL_LEAK = "ROSLEAK_Q_m"
NADPH_OXIDASE = "NOX_YNO1_c"
ROS_REACTION_IDS = (PEROXIDE_EXCHANGE, PEROXIDE_TRANSPORT, MITOCHONDRIAL_LEAK, NADPH_OXIDASE)

EXTRACELLULAR_PEROXIDE = "h2o2_e"
TSA1_ORF = "YML028W"
TSA2_ORF = "YDR453C"
YNO1_ORF = "YGL160W"
THIOREDOXIN_PEROXIDASE = "r_0550"
CYTOSOLIC_CATALASE = "r_0255"
# The Qo site sits inside complex III, so the leak borrows this reaction's GPR verbatim.
COMPLEX_III = "r_0439"

GLUCOSE_EXCHANGE = "r_1714"
OXYGEN_EXCHANGE = "r_1992"
BIOMASS_REACTION = "r_2111"
MAINTENANCE_REACTION = "r_4046"

#: mmol/gDCW/h of H2O2 clearable at glucose 10, oxygen free, CATALASE DELETED. With catalase
#: it is free to a dose of 500 and 1000 still solves, costing 87%: a range, not a property.
REDOX_THROUGHPUT_CEILING = 114.0533

# Every route that consumes cytosolic or nuclear H2O2, i.e. everything a medium dose reaches.
CLEARANCE_REACTIONS = {
    "r_0255": "catalase T, cytosol (CTT1/YGR088W)",
    "r_0483": "glutathione peroxidase, cytosol",
    "r_0550": "thioredoxin peroxidase, cytosol (TSA2, and TSA1 once installed)",
    "r_1037": "thioredoxin peroxidase, nucleus -- structurally blocked, see below",
}

# The three cytosolic routes as gene sets, cheapest first, for :func:`clearance_ladder`.
# r_0256 is out: it is peroxisomal, and no transport takes a cytosolic dose there.
CLEARANCE_ROUTES = (
    ("catalase T", ("YGR088W",), "r_0255"),
    ("glutathione peroxidase",
     ("YNL229C", "YBR244W", "YCL035C", "YDR513W", "YIR037W", "YKL026C"), "r_0483"),
    ("thioredoxin peroxidase", (TSA2_ORF, TSA1_ORF), "r_0550"),
)

# Yeast9 v9.0.2 ids. ecYeastGEM_batch carries all of them as `s_XXXX[c]`, hence _resolve.
_H2O2_C = "s_0837"
_SUPEROXIDE_M, _SUPEROXIDE_C = "s_3813", "s_3931"
_UBIQUINOL_M, _UBIQUINONE_M = "s_1535", "s_1537"
_O2_C, _O2_M = "s_1275", "s_1278"
_PROTON_C, _PROTON_M = "s_0794", "s_0799"
_NADPH_C, _NADP_C = "s_1212", "s_1207"

CITATIONS = {
    PEROXIDE_TRANSPORT: (
        "Bienert & Chaumont 2014, Biochim Biophys Acta 1840:1596-1604; Bienert et al. 2007, "
        "J Biol Chem 282:1183-1192. H2O2 crosses membranes by passive diffusion and through "
        "aquaporins, so an ungated symmetric diffusion is the right shape. Yeast9's own "
        "H2O2 transport r_1839 (cytosol<->nucleus) likewise carries no GPR. NOT EXTENDED to "
        "the mitochondrial or peroxisomal membranes, which the same citation would support: "
        "that would hand a cytosolic dose to the peroxisomal catalase r_0256 and make every "
        "cost here smaller, and it is a change to native compartmentation rather than an "
        "addition, so it is left for a decision rather than taken quietly."
    ),
    PEROXIDE_EXCHANGE: (
        "Boundary reaction, no chemistry to cite. Default bounds (0, 0), shut in BOTH "
        "directions, so installing the module cannot change any existing answer. It was "
        "(0, 1000) and that was wrong: a secretion-only exchange is a free unlimited "
        "peroxide sink, and with it the mutant lacking all nine cytosolic clearance genes "
        "went from dead in the shipped model to 0.887685 /h -- full wild-type growth -- at "
        "glucose 10, dumping 1.69e-4 mmol/gDCW/h of H2O2 into the medium. Only "
        "dose_peroxide opens it now, and it opens it pinned."
    ),
    MITOCHONDRIAL_LEAK: (
        "Murphy 2009, Biochem J 417:1-13; Muller, Liu & Van Remmen 2004, J Biol Chem "
        "279:49064-49073. Univalent reduction of O2 by the complex III Qo-site "
        "ubisemiquinone. The 1 ubiquinol : 2 superoxide ratio is electron conservation "
        "(one electron per superoxide, two per ubiquinol), not a fitted parameter. NO LEAK "
        "FRACTION IS ASSERTED: what fraction of respiratory flux leaks is exactly the "
        "number the ODE layer supplies, and the published estimates disagree by an order of "
        "magnitude, so this module prices a leak and never sets one -- enforced by an upper "
        "bound of ZERO, not merely intended. ASSERTED SIMPLIFICATION: Muller 2004 reports "
        "release to BOTH sides of the inner membrane; all of it is placed in the matrix "
        "here, because that is where SOD2 is and Yeast9 has no intermembrane-space "
        "superoxide. Its protons go to the matrix too, whereas r_0439 puts the Qo-site "
        "protons in the cytosol; that understates the proton-motive force rather than "
        "inventing any, so it is conservative. The GPR is borrowed from r_0439 because the "
        "Qo-site ubisemiquinone needs an assembled complex III, which is an inference from "
        "the mechanism and NOT a finding either citation states. It is also not a "
        "constraint: a GPR limits no flux here, and left open at (0, 1000) this reaction "
        "was an alternative terminal oxidase bypassing COMPLEX IV -- leak, SOD2, then "
        "cytochrome c peroxidase r_0437 -- carrying 43.88 mmol/gDCW/h in a COX1 or COX2 "
        "null, raising that mutant 0.19130 -> 0.37194 /h and doubling its maximum ATP "
        "turnover 20 -> 40, all with complex III present and its GPR intact."
    ),
    NADPH_OXIDASE: (
        "Rinnerthaler et al. 2012, PNAS 109:8658-8663 (Yno1p/Aim14p, YGL160W, is a yeast "
        "NADPH oxidase generating extramitochondrial superoxide); Bedard & Krause 2007, "
        "Physiol Rev 87:245-313 for the NOX one-electron stoichiometry. ASSERTED "
        "SIMPLIFICATION: Yno1p sits in the ER membrane, and the superoxide is placed in the "
        "cytosol because Yeast9 carries no ER superoxide species. Shipped at (0, 0) like "
        "the other two routes: no rate is asserted and the LP may not pick one."
    ),
    TSA1_ORF: (
        "Chae, Chung & Rhee 1994, J Biol Chem 269:27670-27678; Park et al. 2000, J Biol Chem "
        "275:5723-5732; Ogusucu et al. 2007, Free Radic Biol Med 42:326-334 for a "
        "second-order rate constant of order 1e7 /M/s. TSA1 is the dominant cytosolic "
        "peroxiredoxin and was absent from the gene list; NO kinetic parameter enters the "
        "model, the GPR of r_0550 is the entire change."
    ),
}


def _resolve(model: cobra.Model, base_id: str, compartment: str) -> str:
    """Map a bare Yeast9 metabolite id onto this model's convention, or fail by name."""
    for candidate in (base_id, f"{base_id}[{compartment}]", f"{base_id}_{compartment}"):
        if model.metabolites.has_id(candidate):
            return candidate
    raise KeyError(
        f"required metabolite {base_id!r} not found in {model.id!r} "
        f"(tried {base_id}, {base_id}[{compartment}], {base_id}_{compartment})"
    )


def _add_reaction(model, rid, name, stoichiometry, gene_rule="", bounds=(0.0, 1000.0),
                  subsystem="oxidative stress"):
    rxn = cobra.Reaction(rid, name=name, lower_bound=bounds[0], upper_bound=bounds[1])
    rxn.subsystem = subsystem
    model.add_reactions([rxn])
    rxn.add_metabolites({model.metabolites.get_by_id(k): v for k, v in stoichiometry.items()})
    rxn.gene_reaction_rule = gene_rule
    return rxn


def _complex_iii_rule(model: cobra.Model, reaction_id: str = COMPLEX_III) -> str:
    """Complex III's own GPR, read from the model rather than curated here."""
    if not model.reactions.has_id(reaction_id):
        raise KeyError(
            f"{reaction_id} (complex III) is not in {model.id!r}, so the Qo-site leak has "
            "no GPR to borrow and would install as an ungated alternative oxidase"
        )
    rule = model.reactions.get_by_id(reaction_id).gene_reaction_rule
    if not rule.strip():
        raise ValueError(f"{reaction_id} carries no GPR in {model.id!r}")
    return rule


def install_tsa1(model: cobra.Model, reaction_id: str = THIOREDOXIN_PEROXIDASE) -> cobra.Model:
    """Put TSA1 on the cytosolic thioredoxin peroxidase beside TSA2. Modifies ``model``.

    Each TSA2 clause is duplicated with TSA1, because both peroxiredoxins are reduced by
    either thioredoxin. Raises if TSA1 is already on the reaction.
    """
    if not model.reactions.has_id(reaction_id):
        raise KeyError(
            f"{reaction_id} is not in {model.id!r}. GECKO splits it per isoenzyme "
            f"({reaction_id}No1, {reaction_id}No2), and each arm draws on the protein pool, "
            "so adding TSA1 there needs a kcat and a molecular weight. What is published for "
            "Tsa1 is a second-order rate constant, which is not a kcat: refused rather than "
            "converted."
        )
    rxn = model.reactions.get_by_id(reaction_id)
    if TSA1_ORF in rxn.gene_reaction_rule:
        raise ValueError(f"{TSA1_ORF} is already on {reaction_id}: {rxn.gene_reaction_rule}")
    clauses = [c.strip() for c in rxn.gene_reaction_rule.split(" or ")]
    added = [c.replace(TSA2_ORF, TSA1_ORF) for c in clauses if TSA2_ORF in c]
    if not added:
        raise ValueError(
            f"{reaction_id} does not carry TSA2 ({TSA2_ORF}); its rule is "
            f"{rxn.gene_reaction_rule!r}, so there is no paralogue clause to duplicate"
        )
    rxn.gene_reaction_rule = " or ".join(clauses + added)
    model.genes.get_by_id(TSA1_ORF).name = "TSA1"
    return model


def add_ros_module(model: cobra.Model) -> cobra.Model:
    """Return a copy of ``model`` that can be dosed with peroxide and made to leak superoxide.

    The input model is not modified and its maximum growth is unchanged, because every
    addition defaults to a lower bound of zero.

    Raises:
        ValueError: if the module is already installed.
        NotImplementedError: on an enzyme-constrained model.
        KeyError: naming any metabolite this model does not carry.
    """
    if model.reactions.has_id(PEROXIDE_EXCHANGE):
        raise ValueError("the ROS module is already installed on this model")
    # The protein pool is the structural signal that this is a GECKO model, not its name.
    if model.reactions.has_id("prot_pool_exchange"):
        raise NotImplementedError(
            f"{model.id!r} is enzyme-constrained. Every reaction here would carry zero "
            "protein draw, so a GECKO model would run Yno1p and peroxide clearance free of "
            "the pool that is the whole point of it. Installing them honestly needs a kcat "
            "and a molecular weight per enzyme, which this module does not have."
        )

    out = model.copy()
    h2o2_c = _resolve(out, _H2O2_C, "c")
    o2_c, o2_m = _resolve(out, _O2_C, "c"), _resolve(out, _O2_M, "m")
    h_c, h_m = _resolve(out, _PROTON_C, "c"), _resolve(out, _PROTON_M, "m")
    sox_c, sox_m = _resolve(out, _SUPEROXIDE_C, "c"), _resolve(out, _SUPEROXIDE_M, "m")
    qh2, q = _resolve(out, _UBIQUINOL_M, "m"), _resolve(out, _UBIQUINONE_M, "m")
    nadph, nadp = _resolve(out, _NADPH_C, "c"), _resolve(out, _NADP_C, "c")

    out.add_metabolites([
        cobra.Metabolite(EXTRACELLULAR_PEROXIDE, name="hydrogen peroxide",
                         formula="H2O2", charge=0, compartment="e")
    ])
    # Shut BOTH ways by default. Secretion-only is not "no dose", it is a free peroxide sink.
    _add_reaction(out, PEROXIDE_EXCHANGE, "hydrogen peroxide exchange",
                  {EXTRACELLULAR_PEROXIDE: -1}, bounds=(0.0, 0.0), subsystem="exchange")
    _add_reaction(out, PEROXIDE_TRANSPORT, "hydrogen peroxide transport, plasma membrane",
                  {EXTRACELLULAR_PEROXIDE: -1, h2o2_c: 1}, bounds=(-1000.0, 1000.0))
    # Qo-site leak: one electron per superoxide, so one ubiquinol per two. Shut like the
    # exchange, because left open it is an alternative oxidase and S. cerevisiae has none.
    _add_reaction(out, MITOCHONDRIAL_LEAK, "mitochondrial electron leak to superoxide",
                  {qh2: -1, o2_m: -2, q: 1, sox_m: 2, h_m: 2}, bounds=(0.0, 0.0),
                  gene_rule=_complex_iii_rule(out))
    _add_reaction(out, NADPH_OXIDASE, "NADPH oxidase (Yno1p)",
                  {nadph: -1, o2_c: -2, nadp: 1, sox_c: 2, h_c: 1}, bounds=(0.0, 0.0),
                  gene_rule=YNO1_ORF)
    out.genes.get_by_id(YNO1_ORF).name = "YNO1"
    install_tsa1(out)
    return out


def dose_peroxide(model: cobra.Model, dose: float,
                  exchange_id: str = PEROXIDE_EXCHANGE) -> cobra.Model:
    """Force an exact medium peroxide uptake, mmol/gDCW/h. Modifies ``model`` in place.

    Forced rather than merely allowed, because the GEM is a price list: the ODE layer
    supplies the flux and this pins it so the model has to pay for it. Pinning is also what
    makes the dose mandatory work -- merely opening the exchange lets the optimum take up
    whatever suits it, which at glucose 11.1 / oxygen 3.7 is 41.43 mmol/gDCW/h of free
    oxygen. Use inside a ``with model:`` block to scope it.
    """
    if dose < 0:
        raise ValueError(f"dose must be non-negative, got {dose}; it is an uptake magnitude")
    rxn = model.reactions.get_by_id(exchange_id)
    rxn.bounds = (-abs(dose), -abs(dose))
    return model


def dose_superoxide(model: cobra.Model, rate: float,
                    source: str = MITOCHONDRIAL_LEAK) -> cobra.Model:
    """Force an exact superoxide production rate, mmol/gDCW/h. Modifies ``model`` in place.

    Both sources ship at (0, 0) for the same reason the exchange does: left free the LP
    runs them as an alternative oxidase. Use inside a ``with model:`` block to scope it.
    """
    if rate < 0:
        raise ValueError(f"rate must be non-negative, got {rate}; it is a production rate")
    model.reactions.get_by_id(source).bounds = (float(rate), float(rate))
    return model


def _regime(model, glucose_uptake, oxygen_uptake):
    if glucose_uptake is not None:
        model.reactions.get_by_id(GLUCOSE_EXCHANGE).lower_bound = -abs(glucose_uptake)
    if oxygen_uptake is not None:
        model.reactions.get_by_id(OXYGEN_EXCHANGE).lower_bound = -abs(oxygen_uptake)
    return model


def _fluxes(model, reaction_ids, feasible):
    """Flux of each named reaction, or nan when the LP did not solve."""
    return {rid: float(model.reactions.get_by_id(rid).flux) if feasible else float("nan")
            for rid in reaction_ids}


def peroxide_dose_response(
    model: cobra.Model,
    doses,
    glucose_uptake: float | None = 10.0,
    oxygen_uptake: float | None = None,
    knockouts: tuple[str, ...] = (),
    biomass_reaction: str = BIOMASS_REACTION,
) -> pd.DataFrame:
    """Growth against imposed peroxide uptake, with the clearance fluxes that paid for it.

    Args:
        model: A model carrying the ROS module. Not modified.
        doses: Peroxide uptakes to impose, mmol/gDCW/h.
        glucose_uptake: Uptake bound magnitude; ``None`` leaves the model's own.
        oxygen_uptake: Uptake bound magnitude; ``None`` leaves oxygen free.
        knockouts: Gene ORFs to delete first, e.g. the catalases.

    Returns a row per dose with ``growth``, ``feasible``, and the flux through each of
    :data:`CLEARANCE_REACTIONS`. An infeasible dose is reported, not raised.
    """
    if not model.reactions.has_id(PEROXIDE_EXCHANGE):
        raise KeyError("this model has no ROS module; call add_ros_module first")
    rows = []
    for dose in doses:
        with model as m:
            _regime(m, glucose_uptake, oxygen_uptake)
            for orf in knockouts:
                m.genes.get_by_id(orf).knock_out()
            dose_peroxide(m, dose)
            m.objective = biomass_reaction
            value = growth_or_none(m)
            feasible = value is not None
            rows.append({
                "dose": float(dose),
                "growth": value if feasible else float("nan"),
                "feasible": feasible,
                **_fluxes(m, CLEARANCE_REACTIONS, feasible),
            })
    frame = pd.DataFrame(rows)
    reference = frame.loc[frame["dose"] == 0.0, "growth"]
    base = float(reference.iloc[0]) if len(reference) else float("nan")
    frame["growth_cost_fraction"] = 1.0 - frame["growth"] / base
    return frame


def superoxide_leak_response(
    model: cobra.Model,
    leaks,
    source: str = MITOCHONDRIAL_LEAK,
    glucose_uptake: float | None = 10.0,
    oxygen_uptake: float | None = None,
    biomass_reaction: str = BIOMASS_REACTION,
) -> pd.DataFrame:
    """Growth against an imposed superoxide leak, with the SOD flux it forces.

    The superoxide half of the price list. Both SODs were dead reactions before this
    module, so every row here is a number the shipped model could not produce at all.
    """
    if not model.reactions.has_id(source):
        raise KeyError(f"{source} is not in this model; call add_ros_module first")
    rows = []
    for leak in leaks:
        with model as m:
            _regime(m, glucose_uptake, oxygen_uptake)
            m.reactions.get_by_id(source).bounds = (float(leak), float(leak))
            m.objective = biomass_reaction
            value = growth_or_none(m)
            feasible = value is not None
            rows.append({
                "leak": float(leak),
                "growth": value if feasible else float("nan"),
                "feasible": feasible,
                **_fluxes(m, ("r_4190", "r_4270"), feasible),
            })
    frame = pd.DataFrame(rows)
    reference = frame.loc[frame["leak"] == 0.0, "growth"]
    base = float(reference.iloc[0]) if len(reference) else float("nan")
    frame["growth_cost_fraction"] = 1.0 - frame["growth"] / base
    return frame


def clearance_ladder(
    model: cobra.Model,
    dose: float,
    glucose_uptake: float | None = 10.0,
    oxygen_uptake: float | None = None,
    biomass_reaction: str = BIOMASS_REACTION,
) -> pd.DataFrame:
    """Delete the clearance routes cheapest-first and report what the same dose then costs.

    One row per prefix of :data:`CLEARANCE_ROUTES`, starting with none deleted. This is the
    demonstration that a dose is WORK: with every cytosolic route gone the dose is
    infeasible rather than merely expensive, and so is the undosed model, because the
    native oxidases make peroxide that also has to go somewhere.
    """
    rows = []
    for cut in range(len(CLEARANCE_ROUTES) + 1):
        removed = CLEARANCE_ROUTES[:cut]
        orfs = tuple(orf for _label, genes, _rid in removed for orf in genes)
        row = {
            "routes_deleted": cut,
            "deleted": " + ".join(label for label, _g, _r in removed) or "none",
        }
        for key, this_dose in (("undosed_growth", 0.0), ("growth", float(dose))):
            with model as m:
                _regime(m, glucose_uptake, oxygen_uptake)
                for orf in orfs:
                    m.genes.get_by_id(orf).knock_out()
                dose_peroxide(m, this_dose)
                m.objective = biomass_reaction
                value = growth_or_none(m)
                row[key] = value if value is not None else float("nan")
                if key == "growth":
                    row.update(_fluxes(m, [rid for _l, _g, rid in CLEARANCE_ROUTES],
                                       value is not None))
        # A route set whose undosed model cannot grow has no baseline to be a fraction of.
        base = row["undosed_growth"]
        row["growth_cost_fraction"] = (1.0 - row["growth"] / base if base else float("nan"))
        rows.append(row)
    return pd.DataFrame(rows)


@dataclass(frozen=True)
class ClearanceCost:
    """What one clearance route is worth, undosed and dosed.

    ``undosed_cost`` is the standing result that ablating a stress adaptation is free.
    ``dosed_cost`` is what the same deletion costs once a dose has to go somewhere.
    """

    genes: tuple[str, ...]
    dose: float
    glucose_uptake: float | None
    oxygen_uptake: float | None
    undosed_growth: float
    undosed_growth_knocked_out: float
    dosed_growth: float
    dosed_growth_knocked_out: float

    @property
    def undosed_cost(self) -> float:
        return 1.0 - self.undosed_growth_knocked_out / self.undosed_growth

    @property
    def dosed_cost(self) -> float:
        if self.dosed_growth <= 0:
            return float("nan")
        return 1.0 - self.dosed_growth_knocked_out / self.dosed_growth

    def summary(self) -> str:
        regime = (f"glucose {self.glucose_uptake}, oxygen "
                  f"{'free' if self.oxygen_uptake is None else self.oxygen_uptake}")
        return (f"{'+'.join(self.genes)} at {regime}: costs {self.undosed_cost:.3%} of growth "
                f"undosed, {self.dosed_cost:.3%} at a peroxide dose of {self.dose:g} "
                "mmol/gDCW/h")


def clearance_cost(
    model: cobra.Model,
    genes: tuple[str, ...],
    dose: float,
    glucose_uptake: float | None = 10.0,
    oxygen_uptake: float | None = None,
    biomass_reaction: str = BIOMASS_REACTION,
) -> ClearanceCost:
    """Growth cost of deleting a clearance route, with and without a dose to clear.

    This is the whole point of the module in one number: the same deletion that is free in
    an undosed model is not free once peroxide arrives from the medium.
    """

    def growth(d, knock):
        with model as m:
            _regime(m, glucose_uptake, oxygen_uptake)
            if knock:
                for orf in genes:
                    m.genes.get_by_id(orf).knock_out()
            dose_peroxide(m, d)
            m.objective = biomass_reaction
            value = growth_or_none(m)
        return value if value is not None else 0.0

    return ClearanceCost(
        genes=tuple(genes), dose=float(dose),
        glucose_uptake=glucose_uptake, oxygen_uptake=oxygen_uptake,
        undosed_growth=growth(0.0, False), undosed_growth_knocked_out=growth(0.0, True),
        dosed_growth=growth(dose, False), dosed_growth_knocked_out=growth(dose, True),
    )


def energy_from_peroxide(model: cobra.Model, dose: float = 10.0,
                         maintenance_reaction: str = MAINTENANCE_REACTION) -> float:
    """Maximum ATP turnover with peroxide the only thing the cell may take up.

    The classic free-lunch test for a new exchange. A non-zero answer means the module
    built an energy-generating cycle and must not be used.
    """
    with model as m:
        for rxn in m.exchanges:
            rxn.lower_bound = 0.0
        m.reactions.get_by_id(PEROXIDE_EXCHANGE).bounds = (-abs(dose), 0.0)
        maintenance = m.reactions.get_by_id(maintenance_reaction)
        maintenance.bounds = (0.0, 1000.0)
        m.objective = maintenance
        value = growth_or_none(m)
    return 0.0 if value is None else value


def closed_medium_ranges(model: cobra.Model, reaction_ids,
                         maintenance_reaction: str = MAINTENANCE_REACTION) -> pd.DataFrame:
    """FVA with every exchange shut in both directions: the internal-loop free-lunch test.

    Nothing may carry flux when nothing enters the system, so a non-zero range here is a
    thermodynamically infeasible cycle. Maintenance is relaxed first because its non-zero
    lower bound makes a fully closed model infeasible before anything can be measured.
    """
    with model as m:
        for rxn in m.exchanges:
            rxn.bounds = (0.0, 0.0)
        m.reactions.get_by_id(maintenance_reaction).bounds = (0.0, 1000.0)
        m.objective = maintenance_reaction
        return cobra.flux_analysis.flux_variability_analysis(
            m, reaction_list=list(reaction_ids), fraction_of_optimum=0.0,
            processes=FVA_PROCESSES)
