"""Regression contracts for retained legacy kernels and new orchestration."""
import contextlib
import copy
import hashlib
import importlib.util
import io
import tokenize
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch,Mock
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
import workflow_support as w
spec=importlib.util.spec_from_file_location("uq",ROOT/"scripts/04_Uncertainty_Quantification.py")
uq=importlib.util.module_from_spec(spec);spec.loader.exec_module(uq)


class NumericalContract(unittest.TestCase):
    def test_release_sources_match_recorded_hashes(self):
        hashes=json.loads((ROOT/"docs/core_hashes.json").read_text())
        for name,value in hashes.items():
            self.assertEqual(hashlib.sha256((ROOT/name).read_bytes()).hexdigest(),value)

    def test_localized_sources_preserve_original_executable_tokens(self):
        record=json.loads((ROOT/"docs/core_source_localization.json").read_text(encoding="utf8"))
        for name,entry in record["files"].items():
            text=(ROOT/name).read_bytes().decode("utf8")
            try:
                values=[(t.type,t.string) for t in tokenize.generate_tokens(io.StringIO(text).readline) if t.type!=tokenize.COMMENT]
            except tokenize.TokenError:
                values=[(-1,line) for line in text.splitlines(keepends=True) if not line.lstrip().startswith('#')]
            digest=hashlib.sha256(json.dumps(values,ensure_ascii=True,separators=(",",":")).encode()).hexdigest()
            self.assertEqual(digest,entry["non_comment_tokens_sha256"],name)

    def test_entropy_and_ties_are_original_not_normalized(self):
        model=np.array([[1,1],[1,2],[1,3],[1,1],[1,2],[1,3]],dtype=np.uint8)
        with contextlib.redirect_stdout(io.StringIO()): values=w.U.evaluate_nmodel(model,6)
        self.assertEqual(values[0].tolist(),[1,3])
        self.assertAlmostEqual(values[5][0],-2*.0001*np.log(.0001))
        self.assertAlmostEqual(values[5][1],np.log(3))
        self.assertGreater(values[5][1],1.)

    def test_uq_chunking_is_exact(self):
        model=np.random.RandomState(33).randint(1,4,size=(100,19)).astype(np.uint8)
        with tempfile.TemporaryDirectory() as tmp,contextlib.redirect_stdout(io.StringIO()):
            expected=w.U.evaluate_nmodel(model,100)
            actual=uq.legacy_uq(model,Path(tmp),7)
            for name,value in zip(uq.NAMES,expected): np.testing.assert_array_equal(actual[name],value)
            for value in actual.values(): value._mmap.close()

    def test_hierarchical_streaming_keeps_realization_and_voxel_order(self):
        grid=dict(nx=2,ny=2,nz=4,xsiz=10,ysiz=10,zsiz=.5,xmn=5,ymn=5,zmn=.25)
        surfaces=np.array([[.1,.6,1.1,1.6],[1.6,1.1,.6,.1],[.5,.5,1.5,1.5]])
        sis=np.random.RandomState(1).randint(0,2,(3,16))
        with contextlib.redirect_stdout(io.StringIO()),contextlib.redirect_stderr(io.StringIO()):
            expected=w.U.combine_sim_ohnelm(sis,w.U.to_3D(surfaces.ravel(),grid,3),3)
            actual=[w.U.combine_sim_ohnelm(sis[i:i+1],w.U.to_3D(surfaces[i],grid,1),1)[0] for i in range(3)]
        np.testing.assert_array_equal(actual,expected)

    def test_batch_schedule_matches_original_runner(self):
        cfg=w.load_config(ROOT/"configs/smoke.json")
        params=w.sisim_params(cfg,cfg["variograms"]["categorical"])
        calls=[]
        def capture(*args): calls.append(copy.deepcopy(args[-1]))
        with patch.object(w.U,"GSLIB_SISIM",side_effect=capture),patch.object(w.U,"Popen",return_value=Mock()),contextlib.redirect_stdout(io.StringIO()),contextlib.redirect_stderr(io.StringIO()):
            w.U.run_sisim_parallel("sisim.exe","p_","o_","d_","input.gslib",copy.deepcopy(params))
        self.assertEqual(len(calls),20)
        self.assertEqual([p["seed"] for p in calls],[69069+11*i for i in range(20)])
        self.assertTrue(all(p["nsim"]==5 for p in calls))
        self.assertEqual(list(params["grid_params"]),list(w.GRID_KEYS))


if __name__=="__main__": unittest.main()
