<!-- 提交前请确认以下门禁已在本地全部通过；CI 会跑同一组检查 -->

## 本 PR 做了什么 / What this PR does

<!-- 一到三句话说明动机与做法 -->

## 门禁 / Gates（本地已全部通过）

- [ ] `ruff check src tests benchmarks examples testing_suite validation`
- [ ] `ruff format --check src tests benchmarks examples testing_suite validation`
- [ ] `python -m pytest -q`
- [ ] `python testing_suite/scripts/run_all_tests.py`
- [ ] `python examples/quickstart.py --sites 24 --species 60 --trees 3`

## 文档一致性 / Documentation consistency

- [ ] 未改动任何被文档引用的数字，或已同步更新
      （README、docs/、CITATION.cff、.zenodo.json、CHANGELOG）
- [ ] 版本号相关改动已同步 `src/phyloslicer/__init__.py`、
      CITATION.cff、.zenodo.json、README、CHANGELOG

## 类型 / Type

- [ ] bug fix　[ ] feature　[ ] docs　[ ] CI/维护　[ ] 其他
