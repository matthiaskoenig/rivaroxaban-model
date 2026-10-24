"""Annotations of the rivaroxaban models: the model and its species.

The compartments are the shared ones of `pkpdsim.annotations`.
"""

from pkpdsim.annotations import compartments
from sbmlutils.metadata import (
    BQB,
)

__all__ = ["compartments", "model", "species"]

model = [
    # taxonomy
    (BQB.HAS_TAXON, "taxonomy/9606"),  # human
    (BQB.HAS_TAXON, "snomedct/337915000"),  # human
    # modelling approach
    (BQB.HAS_PROPERTY, "mamo/MAMO_0000046"),  # ordinary differential equation model
    # biological process explained by model (GO/NCIT)
    (BQB.HAS_PROPERTY, "ncit/C77995"),  # rivaroxaban
    (BQB.HAS_PROPERTY, "ncit/C180619"),  # Direct Oral Anticoagulant
    (BQB.HAS_PROPERTY, "ncit/C15299"),  # pharmacokinetics
    # model relevance to a particular area
    (BQB.HAS_PROPERTY, "go/GO:0007596"),  # blood coagulation
    # model reference
    (BQB.IS_DESCRIBED_BY, "doi/10.20944/preprints202507.1945.v1"),  # preprint
    (BQB.IS_DESCRIBED_BY, "doi/10.5281/zenodo.15157735"),  # rivaroxaban-model
]

species = {
    "riv": [  # rivaroxaban
        (BQB.IS, "chebi/CHEBI:68579"),
        (BQB.IS, "pubchem.compound/9875401"),
        (BQB.IS, "inchikey/KGFYHTZWPPHNLQ-AWEZNQCLSA-N"),
        (BQB.IS, "ncit/C77995"),
        (BQB.IS, "snomedct/442031002"),
    ],
    # the lumped metabolites of rivaroxaban (M-1, M-2, M-4, M-7 and others,
    # Weinz2009), which have no identifiers of their own: metabolites
    "rx": [
        (BQB.IS_VERSION_OF, "ncit/C1513"),  # metabolite
    ],
}
