"""Guard the scientific default, historical CLI and frozen-input identity."""

import pytest

from ensemble import benchmark, replay_selected


def test_default_and_explicit_legacy_routes(monkeypatch):
    invoked = []
    monkeypatch.setattr(replay_selected, 'main', lambda args: invoked.append(('full5', args)))
    monkeypatch.setattr(benchmark, 'legacy_main', lambda args: invoked.append(('legacy', args)))
    benchmark.main(['--output-root', 'new'])
    benchmark.main(['--legacy-four', '--help'])
    benchmark.main(['--artifact-root=old', '--output-root', 'historical'])
    assert invoked == [('full5', ['--output-root', 'new']), ('legacy', ['--help']),
                       ('legacy', ['--artifact-root=old', '--output-root', 'historical'])]


def test_unregistered_input_manifest_does_not_publish_a_result(tmp_path):
    inputs = tmp_path / 'inputs'
    inputs.mkdir()
    (inputs / 'manifest.json').write_text('{}')
    output = tmp_path / 'scores'
    with pytest.raises(ValueError, match='registered'):
        replay_selected.main(['--input-root', str(inputs), '--output-root', str(output)])
    assert not output.exists()


def test_existing_output_is_preserved(tmp_path):
    marker = tmp_path / 'keep.txt'
    marker.write_text('previous run')
    with pytest.raises(SystemExit, match='2'):
        replay_selected.main(['--output-root', str(tmp_path)])
    assert marker.read_text() == 'previous run'


def test_historical_build_cannot_replace_selected_data():
    selected = replay_selected.ROOT / "ensemble/results/selected"
    with pytest.raises(SystemExit, match='2'):
        benchmark.main(['--artifact-root', 'historical-inputs', '--output-root', str(selected)])
