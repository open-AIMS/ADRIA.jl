"""Small deterministic checks for the 1991 historical surface builder."""

from __future__ import annotations

import unittest

import numpy as np
import pandas as pd

from build_1991_initial_surfaces import crosswalk, haversine_km, idw


class HistoricalSurfaceTests(unittest.TestCase):
    def test_effort_weighting_and_missing_coverage(self) -> None:
        distances = np.array([[1.0, 2.0, 3.0], [1.0, 2.0, 120.0]])
        values = np.array([0.0, 0.2, 0.4])
        precision = np.array([10.0, 20.0, 10.0])
        estimate, donors, _ = idw(distances, values, precision)
        expected = (10 * 0 + 20 / 4 * 0.2 + 10 / 9 * 0.4) / (10 + 20 / 4 + 10 / 9)
        self.assertAlmostEqual(estimate[0], expected)
        self.assertEqual(donors.tolist(), [3, 2])
        self.assertTrue(np.isnan(estimate[1]))

    def test_haversine_zero(self) -> None:
        self.assertAlmostEqual(haversine_km(np.array([145.0]), np.array([-14.0]),
                                             np.array([145.0]), np.array([-14.0]))[0, 0], 0.0)

    def test_ambiguous_name_is_not_guessed(self) -> None:
        reefs = pd.DataFrame({"LABEL_ID": ["14-001", "14-002"],
                              "reefName": ["Same Reef (14-001)", "Same Reef (14-002)"],
                              "Shape_Area": [1.0, 1.0],
                              "xCentroid": [145.0, 145.1],
                              "yCentroid": [-14.0, -14.1]})
        result = crosswalk(reefs, {"Same Reef"})
        self.assertEqual(result.status.iloc[0], "ambiguous")


if __name__ == "__main__":
    unittest.main()
