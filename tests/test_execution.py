from pathlib import Path

from abag_rbfe.execution import CommandRunner


def test_command_runner_streams_success_output_to_log(tmp_path: Path) -> None:
    runner = CommandRunner(execute=True)
    script_path = tmp_path / "artifacts" / "commands" / "demo.sh"

    outcome = runner.run_script(
        script_path,
        commands=[
            "echo first-line",
            "echo second-line",
        ],
        workdir=tmp_path,
    )

    log_path = script_path.with_suffix(".log")
    assert outcome.state == "completed"
    assert str(log_path) in outcome.message
    assert log_path.read_text(encoding="utf-8").splitlines() == ["first-line", "second-line"]


def test_command_runner_writes_failure_output_to_log(tmp_path: Path) -> None:
    runner = CommandRunner(execute=True)
    script_path = tmp_path / "artifacts" / "commands" / "fail.sh"

    outcome = runner.run_script(
        script_path,
        commands=[
            "echo before-failure",
            "bash -c 'echo fatal-message >&2; exit 3'",
        ],
        workdir=tmp_path,
    )

    log_path = script_path.with_suffix(".log")
    log_text = log_path.read_text(encoding="utf-8")
    assert outcome.state == "failed"
    assert str(log_path) in outcome.message
    assert "before-failure" in log_text
    assert "fatal-message" in log_text
    assert "exit code 3" in outcome.message


def test_command_runner_reports_signal_termination(tmp_path: Path) -> None:
    runner = CommandRunner(execute=True)
    script_path = tmp_path / "artifacts" / "commands" / "term.sh"

    outcome = runner.run_script(
        script_path,
        commands=[
            "echo before-term",
            "bash -c 'kill -TERM $$'",
        ],
        workdir=tmp_path,
    )

    log_path = script_path.with_suffix(".log")
    assert outcome.state == "failed"
    assert str(log_path) in outcome.message
    assert "signal 15" in outcome.message
    assert "before-term" in log_path.read_text(encoding="utf-8")
