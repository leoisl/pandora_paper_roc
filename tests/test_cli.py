from evaluate.cli import *
import pytest
import uuid


class TestParseArgs:
    def test_noArgsGivenRaisesError(self):
        args = []

        with pytest.raises(SystemExit):
            parse_args(args)

    def test_invalidOutputGivenRaisesError(self):
        args = "-1 foo -2 bar -v ni -r blah -o /fakedir/out.tsv".split()

        with pytest.raises(NotADirectoryError):
            parse_args(args)

    def test_invalidLogLevelRaisesError(self):
        log_lvl = 6
        args = f"-1 foo -2 bar -v ni -r blah --log_level {log_lvl}".split()

        with pytest.raises(SystemExit):
            parse_args(args)

    def test_tmpDirCreated(self):
        temp = Path(f"/tmp/{str(uuid.uuid4())}")
        args = f"-1 foo -2 bar -v ni -r blah --temp {str(temp)}".split()

        assert not temp.exists()

        parse_args(args)

        assert temp.exists()

    def test_outputFileOpenedForWriting(self):
        output = Path(f"/tmp/{str(uuid.uuid4())}")
        args = f"-1 foo -2 bar -v ni -r blah --output {str(output)}".split()

        assert not output.exists()

        actual = parse_args(args)

        assert actual.output.name == str(output)
        assert not actual.output.closed

        actual.output.close()
        output.unlink()
