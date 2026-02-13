"""
Test database URLs to ensure they are reachable and valid.

This test checks that all URLs in downloads.json are accessible without
downloading the full files (which can be very large).
"""

import json
import os

import pytest
import requests


@pytest.fixture(scope="module")
def database_urls():
    """Load database URLs from downloads.json."""
    downloads_json = os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(__file__))),
        "funannotate2",
        "downloads.json",
    )
    with open(downloads_json, "r") as f:
        data = json.load(f)
    return data["downloads"]


def check_url_reachable(url, timeout=10):
    """
    Check if a URL is reachable without downloading the full content.

    Args:
        url: The URL to check
        timeout: Request timeout in seconds

    Returns:
        tuple: (is_reachable: bool, status_code: int, error_message: str or None)
    """
    try:
        # Try HEAD request first (faster, doesn't download content)
        response = requests.head(url, timeout=timeout, allow_redirects=True)

        # Some servers don't support HEAD, so if we get 405 or 501, try GET with range
        if response.status_code in [405, 501]:
            # Try GET request with small range to avoid downloading full file
            headers = {"Range": "bytes=0-1023"}  # Only get first 1KB
            response = requests.get(
                url, headers=headers, timeout=timeout, allow_redirects=True, stream=True
            )
            # Close the connection immediately
            response.close()

        if response.status_code in [200, 206, 302, 301]:  # 206 = Partial Content
            return True, response.status_code, None
        else:
            return False, response.status_code, f"HTTP {response.status_code}"

    except requests.exceptions.Timeout:
        return False, None, f"Timeout after {timeout}s"
    except requests.exceptions.ConnectionError as e:
        return False, None, f"Connection error: {str(e)}"
    except requests.exceptions.RequestException as e:
        return False, None, f"Request failed: {str(e)}"
    except Exception as e:
        return False, None, f"Unexpected error: {str(e)}"


class TestDatabaseURLs:
    """Test that database URLs are reachable."""

    @pytest.mark.parametrize(
        "db_name",
        [
            "uniprot",
            "uniprot-release",
            "merops",
            "dbCAN",
            "dbCAN-tsv",
            "dbCAN-log",
            "pfam",
            "pfam-tsv",
            "pfam-log",
            "go",
            "mibig",
            "interpro",
            "interpro-tsv",
            "gene2product",
            "mito",
            "mito-release",
        ],
    )
    def test_database_url_reachable(self, database_urls, db_name):
        """Test that each database URL is reachable."""
        url = database_urls.get(db_name)
        assert url is not None, f"URL for {db_name} not found in downloads.json"

        is_reachable, status_code, error_msg = check_url_reachable(url, timeout=15)

        assert is_reachable, (
            f"URL for {db_name} is not reachable:\n"
            f"  URL: {url}\n"
            f"  Status: {status_code}\n"
            f"  Error: {error_msg}"
        )

    def test_all_urls_present(self, database_urls):
        """Test that all expected database URLs are present in downloads.json."""
        expected_dbs = [
            "uniprot",
            "uniprot-release",
            "merops",
            "dbCAN",
            "dbCAN-tsv",
            "dbCAN-log",
            "pfam",
            "pfam-tsv",
            "pfam-log",
            "go",
            "mibig",
            "interpro",
            "interpro-tsv",
            "gene2product",
            "mito",
            "mito-release",
        ]

        for db_name in expected_dbs:
            assert (
                db_name in database_urls
            ), f"Expected database '{db_name}' not found in downloads.json"

    def test_urls_are_valid_format(self, database_urls):
        """Test that all URLs are in valid format."""
        for db_name, url in database_urls.items():
            assert isinstance(url, str), f"URL for {db_name} is not a string"
            assert url.startswith(("http://", "https://", "ftp://")), (
                f"URL for {db_name} does not start with http://, https://, or ftp://: {url}"
            )

