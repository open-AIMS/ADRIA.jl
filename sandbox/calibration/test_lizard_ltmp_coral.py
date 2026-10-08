"""Coverage and alias checks for the corrected LTMP coral diagnostic."""

import unittest

from plot_lizard_ltmp_coral import ROOT, REEF_ALIASES, load_ltmp_coral


class LTMPCoralMappingTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.coral = load_ltmp_coral(
            ROOT / "sandbox/data/reef_manta.csv",
            ROOT / "sandbox/data/reef_photo_transect.csv",
        )

    def test_manta_covers_all_four_aliases(self):
        counts = self.coral[self.coral.method == "manta"].groupby("model_reef_name").size()
        self.assertEqual(set(counts.index), set(REEF_ALIASES))
        self.assertEqual(counts.to_dict(), {
            "Lizard Island Reef": 32,
            "MacGillivray Reef": 33,
            "North Direction Reef": 27,
            "Eyrie Reef": 11,
        })

    def test_photo_uses_group_level_and_excludes_eyrie(self):
        counts = self.coral[self.coral.method == "photo_transect"].groupby("model_reef_name").size()
        self.assertEqual(counts.to_dict(), {
            "Lizard Island Reef": 23,
            "MacGillivray Reef": 24,
            "North Direction Reef": 23,
        })

    def test_2015_reference_values_and_unique_keys(self):
        subset = self.coral[self.coral.year == 2015]
        self.assertEqual(len(subset), 7)
        self.assertFalse(self.coral.duplicated(["model_reef_name", "method", "year"]).any())
        expected = {
            ("Lizard Island Reef", "manta"): 0.0788522250,
            ("MacGillivray Reef", "manta"): 0.0842473510,
            ("North Direction Reef", "manta"): 0.0676574233,
            ("Eyrie Reef", "manta"): 0.1046277473,
            ("Lizard Island Reef", "photo_transect"): 0.0854283254,
            ("MacGillivray Reef", "photo_transect"): 0.1181926954,
            ("North Direction Reef", "photo_transect"): 0.0883053341,
        }
        for row in subset.itertuples():
            self.assertAlmostEqual(row.median, expected[(row.model_reef_name, row.method)], places=8)


if __name__ == "__main__":
    unittest.main()
