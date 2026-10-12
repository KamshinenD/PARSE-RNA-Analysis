"""Verify the direct-input contract and timing boundary."""
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from types import SimpleNamespace

import benchmark_direct as direct


class DirectTests(unittest.TestCase):
    def test_native_commands_have_no_preprocessing_or_wrapper(self):
        config = dict(python="/venv/bin/python", parse="/bin/parse", dssr="/bin/dssr",
                      rnaview="/bin/rnaview", bpnet="/bin/bpnet")
        work = Path("/tmp/run")
        cif = work / "TEST.cif"
        self.assertEqual(direct.command_for("dssr", cif, work, config),
                         ["/bin/dssr", "-i=/tmp/run/TEST.cif", "--pair-only"])
        self.assertEqual(direct.command_for("parse", cif, work, config),
                         ["/bin/parse", str(cif), "--no-download", "--cache-dir", str(work)])
        self.assertEqual(direct.command_for("barnaba", cif, work, config),
                         ["/venv/bin/barnaba", "ANNOTATE", "--pdb", str(cif), "-o", "annotation"])
        for tool in direct.common.TOOLS:
            self.assertNotIn("run_one.py", " ".join(direct.command_for(tool, cif, work, config)))

    def test_original_bytes_and_validation_outside_stopwatch(self):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            (folder / "logs").mkdir()
            source = folder / "TEST.cif"
            original = b"data_original\n# Keep every model and chain exactly as supplied\n"
            source.write_bytes(original)
            events = []
            def launch(command, **kwargs):
                cif = Path(command[1])
                self.assertTrue(cif.is_symlink())
                self.assertEqual(cif.resolve(), source.resolve())
                self.assertEqual(cif.read_bytes(), original)
                events.append("launch")
                def communicate(**kwargs):
                    events.append("exit")
                    return b'{"n_pairs":0,"pairs":[]}', b""
                return SimpleNamespace(returncode=0, communicate=communicate)
            def clock():
                events.append("clock")
                return 10.0 if len(events) == 1 else 12.5
            def validate(*args):
                events.append("validate")
                return []
            with patch.object(direct.subprocess, "Popen", side_effect=launch), \
                 patch.object(direct.time, "perf_counter", side_effect=clock), \
                 patch.object(direct, "validate_output", side_effect=validate):
                result = direct.run_trial("parse", dict(pdb="TEST", nucleotides=1, source_path=str(source)),
                    0, SimpleNamespace(timeout=10), {"parse":"parse", "python":"/venv/bin/python"}, folder, folder, "originalhash")
            self.assertEqual(events, ["clock", "launch", "exit", "clock", "validate"])
            self.assertEqual(result["wall_seconds"], 2.5)
            self.assertEqual(result["status"], "ok")
            self.assertEqual(source.read_bytes(), original)

    def test_recorded_command_and_config_have_no_machine_paths(self):
        config = dict(python="/home/me/venv/bin/python", parse="/home/me/build/parse",
                      dssr="/opt/x/x3dna-dssr", rnaview="/home/me/rv/bin/rnaview",
                      bpnet="/home/me/bp/bpnet", packages={"numpy": "2.2.6"})
        work = Path("/tmp/rnad-abc")
        cif = work / "1ABC.cif"
        for tool in direct.common.TOOLS:
            shown = direct.public_command(direct.command_for(tool, cif, work, config), work, config)
            self.assertFalse(any(part.startswith("/") for part in shown), (tool, shown))
        self.assertEqual(direct.public_command(direct.command_for("dssr", cif, work, config), work, config),
                         ["x3dna-dssr", "-i=WORK/1ABC.cif", "--pair-only"])
        public = direct.public_config(config)
        self.assertFalse(any(isinstance(v, str) and v.startswith("/") for v in public.values()))
        self.assertEqual(public["packages"], {"numpy": "2.2.6"})

    def test_truncated_dssr_output_rejected(self):
        with self.assertRaises(ValueError):
            direct.validate_output("dssr", Path("x.cif"), Path("."), "", "total number of base pairs: 2")
        self.assertEqual(direct.validate_output("dssr", Path("x.cif"), Path("."), "",
                                              "total number of base pairs: 0"), [])


if __name__ == "__main__":
    unittest.main()
