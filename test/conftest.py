import os
import sys
from types import SimpleNamespace

import pytest


this_dir = os.path.dirname(__file__)
DATA = os.path.join(this_dir, "data")


@pytest.fixture(scope="session")
def dec():
    """
    The decifer modules, imported with a command line pointing at test/data: decifer.parse_args
    and decifer.process_input parse sys.argv and read the input when first imported.
    """
    argv = sys.argv
    sys.argv = ["decifer", os.path.join(DATA, "decifer_input.tsv"), "-p", os.path.join(DATA, "purity.tsv")]
    try:
        import decifer.parse_args
        import decifer.process_input
        import decifer.config
        import decifer.mutation
        import decifer.fileio
        import decifer.new_coordinate_ascent
        import decifer.fastmodel
    finally:
        sys.argv = argv
    state_trees = decifer.fileio.read_in_state_trees(decifer.parse_args.args["statetrees"])
    return SimpleNamespace(
        args=decifer.parse_args.args,
        purity=decifer.process_input.PURITY,
        df=decifer.process_input.MUTATION_DF,
        state_trees=state_trees,
        config=decifer.config,
        mutation=decifer.mutation,
        nca=decifer.new_coordinate_ascent,
        fastmodel=decifer.fastmodel,
    )
