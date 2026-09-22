"""From a reaction the agent chose to the gene edits a laboratory would actually make.

The agent picks reactions, because that is what flux-based evidence is about: FSEOF ranks
reactions, OptKnock deletes reactions, and a distance to the product is a distance through
reactions. But nobody edits a reaction. What gets built is a gene deletion or a weakened
promoter, and the map from one to the other is the gene-protein-reaction rule — which is not
a list of names, and is not a bijection.

Three ways the naive reading goes wrong, all of them common:

**Isozymes.** ``ACKr`` in ``e_coli_core`` is ``b2296 or b3115 or b1849``. Deleting *ackA*, the
textbook acetate-branch gene, leaves two other routes and changes nothing at all — measured,
succinate stays at 0.0000 and growth at 0.2117. To block the reaction you must delete all
three. On this model 32 of the 69 gene-associated reactions are like this.

**Complexes.** ``THD2`` is ``b1602 and b1603``. Here the opposite holds: deleting *either*
subunit is enough, so naming both overstates the work. Ten reactions.

**Shared genes.** A gene serves whatever reactions it serves. Deleting all of ``ACALDt``'s
genes also blocks ``CO2t`` and ``O2t``; ``SUCCt2_2`` shares its transporter with ``FUMt2_2``
and ``MALt2_2``. Nineteen reactions on this model take others down with them, and a bound set
on one reaction shows none of it.

So a move is *resolved* here before it is applied: the minimal set of genes whose loss stops
the chosen reaction, and every other reaction that set stops too. The engine then applies the
whole consequence, which is what makes CMM's viability rule a rule about strains that can be
built rather than about bound edits that only exist in silico.
"""

from __future__ import annotations

import ast
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass

from cobra import Model, Reaction


@dataclass(frozen=True)
class GeneEdit:
    """What editing one reaction actually means: which genes, and what else it hits.

    ``genes`` is a minimal set — remove any one of them and the target reaction comes back —
    so it is the smallest piece of laboratory work that achieves the move, not every gene the
    reaction mentions.

    ``blocks`` is every reaction the edit stops, the target included. A caller that applies
    bounds to ``target`` alone is modelling a strain that cannot be built.
    """

    target: str
    genes: tuple[str, ...]
    blocks: tuple[str, ...]
    #: True when the rule is a pure OR of single genes, so each one alone is useless. Worth
    #: saying out loud: it is the case people get wrong, and it changes how much work the move
    #: is by a factor of three on ``ACKr``.
    isozymes: bool = False
    #: True when some OR-term needed more than one gene, so alternatives existed and the
    #: cheapest was taken.
    chose_among_alternatives: bool = False

    @property
    def side_effects(self) -> tuple[str, ...]:
        """Reactions stopped that the agent did not choose."""

        return tuple(r for r in self.blocks if r != self.target)

    def describe(self, names: Mapping[str, str] | None = None) -> str:
        """A phrase a wet-lab reader can act on, naming genes rather than loci where it can."""

        names = names or {}

        def label(gene: str) -> str:
            return names.get(gene) or gene

        if not self.genes:
            return "no gene association; this is a bound edit, not a gene edit"
        listed = ", ".join(label(g) for g in self.genes)
        if len(self.genes) == 1:
            phrase = f"gene {listed}"
        elif self.isozymes:
            phrase = f"all of {listed} together (isozymes: any one alone does nothing)"
        else:
            phrase = f"{listed} together"
        if self.side_effects:
            phrase += f" — which also stops {', '.join(self.side_effects)}"
        return phrase


def _or_terms(node: ast.expr | None) -> list[frozenset[str]]:
    """The GPR body as a disjunction of gene sets: ``(a and b) or c`` -> ``[{a,b}, {c}]``.

    Blocking the reaction means breaking every term, so the terms are what a gene set has to
    hit. Nested rules are flattened by distributing ``and`` over ``or``; GPRs are small enough
    that the expansion is never a problem in practice.
    """

    if node is None:
        return []
    if isinstance(node, ast.BoolOp) and isinstance(node.op, ast.Or):
        terms: list[frozenset[str]] = []
        for value in node.values:
            terms.extend(_or_terms(value))
        return terms
    if isinstance(node, ast.BoolOp) and isinstance(node.op, ast.And):
        combined: list[frozenset[str]] = [frozenset()]
        for value in node.values:
            expanded = [base | term for term in _or_terms(value) for base in combined]
            combined = expanded or combined
        return combined
    if isinstance(node, ast.Name):
        return [frozenset({node.id})]
    return []  # pragma: no cover - cobra does not produce other node types here


def _reachable(model: Model, genes: Iterable[str]) -> list[Reaction]:
    """Every reaction any of these genes takes part in.

    The collateral of a gene edit can only be here, so the GPR of everything else never has to
    be evaluated. On a genome-scale model that is the difference between a handful of checks
    and one per reaction, on a path that runs twice per candidate per screen.
    """

    touched: dict[str, Reaction] = {}
    for gene_id in genes:
        try:
            gene = model.genes.get_by_id(gene_id)
        except KeyError:  # pragma: no cover - the ids come from the model itself
            continue
        for reaction in gene.reactions:
            touched[reaction.id] = reaction
    return list(touched.values())


def resolve_gene_edit(model: Model, reaction_id: str) -> GeneEdit:
    """The smallest gene deletion that stops this reaction, and everything else it stops.

    Minimal in cardinality: one gene is needed from every OR-term and no more, so the set is
    as small as it can be. Among the sets of that size the one chosen does the least collateral
    damage — greedily, by preferring at each step the gene that breaks the most remaining terms
    and, among equals, the one serving the fewest other reactions. Ties then break on the gene
    id, so a run is reproducible.

    The result is *checked* against cobra's own GPR evaluation rather than trusted from the
    parse, because a rule this function read wrongly would otherwise become a design nobody
    could build.
    """

    reaction = model.reactions.get_by_id(reaction_id)
    gpr = reaction.gpr
    terms = _or_terms(getattr(gpr, "body", None))
    if not terms or not reaction.genes:
        return GeneEdit(target=reaction_id, genes=(), blocks=(reaction_id,))

    # How many other reactions each gene would take down with it, used only to break ties.
    cost = {
        gene.id: len([r for r in gene.reactions if r.id != reaction_id])
        for gene in reaction.genes
    }

    chosen: list[str] = []
    remaining = list(terms)
    alternatives = any(len(term) > 1 for term in terms)
    while remaining:
        candidates = sorted(
            {gene for term in remaining for gene in term},
            key=lambda gene: (
                -sum(1 for term in remaining if gene in term),
                cost.get(gene, 0),
                gene,
            ),
        )
        pick = candidates[0]
        chosen.append(pick)
        remaining = [term for term in remaining if pick not in term]

    # Drop anything the greedy pass made redundant, cheapest-to-lose first.
    minimal = list(chosen)
    for gene in sorted(chosen, key=lambda g: (-cost.get(g, 0), g)):
        trial = [g for g in minimal if g != gene]
        if trial and gpr.eval(trial) is False:
            minimal = trial

    genes = tuple(sorted(minimal))
    if not genes or gpr.eval(list(genes)) is not False:  # pragma: no cover - defensive
        # The parse disagreed with cobra. Fall back to deleting everything the reaction
        # names, which is always sufficient, and say nothing clever about minimality.
        genes = tuple(sorted(g.id for g in reaction.genes))

    blocked = {reaction_id}
    for other in _reachable(model, genes):
        if other.gpr.eval(list(genes)) is False:
            blocked.add(other.id)
    return GeneEdit(
        target=reaction_id,
        genes=genes,
        blocks=tuple(sorted(blocked)),
        isozymes=len(terms) > 1 and all(len(term) == 1 for term in terms),
        chose_among_alternatives=alternatives,
    )


def gene_names(model: Model, genes: Sequence[str]) -> dict[str, str]:
    """Locus tag -> the readable name, where the model carries one.

    A brief that says "delete ldhA" and a board that says ``b1380`` are not connectable by
    anything, which was the whole point of putting genes on the board.
    """

    names: dict[str, str] = {}
    for gene_id in genes:
        try:
            gene = model.genes.get_by_id(gene_id)
        except KeyError:  # pragma: no cover - ids come from the model
            continue
        name = str(getattr(gene, "name", "") or "").strip()
        if name and name != gene_id:
            names[gene_id] = name
    return names


__all__ = ["GeneEdit", "gene_names", "resolve_gene_edit"]
