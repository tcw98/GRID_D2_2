"""Catch inconsistent CV tables or accidental replacement of the selected input."""
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from verify_final_validation import verify


class FinalValidationEvidence(unittest.TestCase):
    def test_packaged_cv_and_spatial_evidence(self):
        result = verify(ROOT)
        self.assertEqual(result['native_realizations'], 500)
        self.assertEqual(result['binary_N'], 22313)
        self.assertEqual(result['continuous_N'], 664)


if __name__ == '__main__':
    unittest.main()
