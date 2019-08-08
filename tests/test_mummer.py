import pytest
import uuid
from pathlib import Path
from evaluate.mummer import Nucmer, DeltaFilter

cft = Path("tests/test_cases/CFT073.ref.fa")
h13 = Path("tests/test_cases/H131800734.ref.pilon.fa")
delta = Path("tests/test_cases/out.delta")


class TestNucmer:
    def test_init_withFakeFileRaisesFileNotFoundError(self):
        reference = Path("foo")
        query = Path("bar")
        with pytest.raises(FileNotFoundError):
            nucmer = Nucmer(reference, query)

    def test_generateCommand_defaultArgsHasNoExtraParamsAndHasOutAsPrefix(self):
        reference = cft
        query = h13
        nucmer = Nucmer(reference, query)

        actual = nucmer.generate_command()
        expected = ["nucmer", "--prefix", "out", str(reference), str(query)]

        assert actual == expected

    def test_generateCommand_extraParamsIncluded(self):
        reference = cft
        query = h13
        extra_params = "--maxmatch --forward"
        nucmer = Nucmer(reference, query, extra_params=extra_params)

        actual = nucmer.generate_command()
        expected = [
            "nucmer",
            extra_params,
            "--prefix",
            "out",
            str(reference),
            str(query),
        ]

        assert actual == expected

    @pytest.mark.skip(reason="Takes a while to run so don't want to run all the time.")
    def test_run_twoTestCases_returnsOkExitCode(self):
        reference = cft
        query = h13
        prefix = f"/tmp/{str(uuid.uuid4())}"
        extra_params = "--maxmatch"
        nucmer = Nucmer(reference, query, prefix=prefix, extra_params=extra_params)
        result = nucmer.run()

        actual = result.returncode
        expected = 0

        assert actual == expected

        Path(prefix + ".delta").unlink()


class TestDeltaFilter:
    def test_init_withFakeFileRaisesFileNotFoundError(self):
        deltafile = Path("foo")

        with pytest.raises(FileNotFoundError):
            deltafilter = DeltaFilter(deltafile)

    def test_generateCommand_defaultArgsHasNoExtraParamsAndHasOutAsPrefix(self):
        deltafile = delta
        deltafilter = DeltaFilter(deltafile)

        actual = deltafilter.generate_command()
        expected = ["delta-filter", str(deltafile)]

        assert actual == expected

    def test_generateCommand_extraParamsIncluded(self):
        deltafile = delta
        extra_params = "-1"
        deltafilter = DeltaFilter(deltafile, extra_params=extra_params)

        actual = deltafilter.generate_command()
        expected = ["delta-filter", extra_params, str(deltafile)]

        assert actual == expected

    def test_run_realDeltaFileReturnsOkExitCode(self):
        deltafile = delta
        extra_params = "-1"
        deltafilter = DeltaFilter(deltafile, extra_params=extra_params)
        result = deltafilter.run()

        actual = result.returncode
        expected = 0

        assert actual == expected

        actual = result.stdout.decode().split("\n")[2:4]
        expected = [
            ">1 1_pilon_pilon_pilon_pilon_pilon_pilon 5155066 4725092",
            "1238 2091 273232 274082 62 62 0",
        ]

        assert actual == expected
