"""Tests for ancestry estimation module."""

from genomeinsight.ancestry.reference_data import AIMDatabase


class TestAIMDatabase:
    """Tests for AIM reference data loading."""

    def test_load_returns_database(self):
        db = AIMDatabase.load()
        assert db is not None

    def test_has_populations(self):
        db = AIMDatabase.load()
        assert len(db.populations) == 14

    def test_has_aims(self):
        db = AIMDatabase.load()
        assert len(db.aims) >= 100

    def test_population_has_region(self):
        db = AIMDatabase.load()
        for pop_name, pop_info in db.populations.items():
            assert "region" in pop_info, f"{pop_name} missing region"

    def test_aim_has_frequencies_for_all_populations(self):
        db = AIMDatabase.load()
        pop_names = set(db.populations.keys())
        for rsid, aim in db.aims.items():
            aim_pops = set(aim["frequencies"].keys())
            assert aim_pops == pop_names, f"{rsid} missing populations: {pop_names - aim_pops}"

    def test_frequencies_are_valid(self):
        db = AIMDatabase.load()
        for rsid, aim in db.aims.items():
            for pop, freq in aim["frequencies"].items():
                assert 0.0 <= freq <= 1.0, f"{rsid} {pop}: freq {freq} out of range"

    def test_get_aim_rsids(self):
        db = AIMDatabase.load()
        rsids = db.get_aim_rsids()
        assert isinstance(rsids, set)
        assert len(rsids) == len(db.aims)

    def test_get_population_names(self):
        db = AIMDatabase.load()
        names = db.get_population_names()
        assert "Northern European" in names
        assert "West African" in names
