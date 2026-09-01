"""
Unit tests for functional-annotation key resolution in utilities.py.
"""

from funannotate2.utilities import resolve_annotation_key


class TestResolveAnnotationKey:
    """Regression coverage for issue #96.

    annotate() merged functional annotations with exact-key checks
    (`if "ec_number" in fa` / `if "go_terms" in fa`), but eggNOG/addons emit the
    EC number under `EC_number` and InterProScan emits GO terms under `go_term`
    (singular), so those annotations were silently dropped. resolve_annotation_key
    tolerates the known case/plural variants.
    """

    def test_canonical_ec_number(self):
        assert (
            resolve_annotation_key({"ec_number": [["1.1.1.1"]]}, "ec_number")
            == "ec_number"
        )

    def test_capitalized_ec_number_variant(self):
        # eggNOG / funannotate2-addons emit "EC_number"; dropped before the fix.
        assert (
            resolve_annotation_key({"EC_number": [["1.1.1.1"]]}, "ec_number")
            == "EC_number"
        )

    def test_canonical_go_terms(self):
        assert (
            resolve_annotation_key({"go_terms": [["GO:0008150"]]}, "go_terms")
            == "go_terms"
        )

    def test_singular_go_term_variant(self):
        # InterProScan emits "go_term" (singular); dropped before the fix.
        assert (
            resolve_annotation_key({"go_term": [["GO:0008150"]]}, "go_terms")
            == "go_term"
        )

    def test_missing_key_returns_none(self):
        assert resolve_annotation_key({"product": ["x"]}, "ec_number") is None
        assert resolve_annotation_key({}, "go_terms") is None

    def test_canonical_preferred_when_both_present(self):
        # Canonical wins so canonical-key inputs behave exactly as before #96.
        assert (
            resolve_annotation_key({"ec_number": 1, "EC_number": 2}, "ec_number")
            == "ec_number"
        )

    def test_unknown_canonical_key_matches_itself_only(self):
        # A key with no known variants resolves only to itself.
        assert resolve_annotation_key({"note": ["x"]}, "note") == "note"
        assert resolve_annotation_key({"other": 1}, "note") is None

    def test_populated_variant_wins_over_empty_canonical(self):
        # #96 intent is to never re-drop data: an empty canonical key must not
        # shadow a populated variant when both are present.
        assert (
            resolve_annotation_key(
                {"ec_number": [], "EC_number": [["1.1.1.1"]]}, "ec_number"
            )
            == "EC_number"
        )

    def test_present_but_none_value_preserved(self):
        # A present key with a falsy value still resolves (the old code used
        # key-presence, not truthiness); the annotate() call site then assigns it,
        # preserving the pre-fix behavior for that edge.
        assert resolve_annotation_key({"ec_number": None}, "ec_number") == "ec_number"
