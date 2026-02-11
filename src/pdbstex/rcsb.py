from __future__ import annotations

import json
import time
from dataclasses import dataclass
from typing import Any, Optional

import requests

from .utils import WarningCollector, sleep_backoff, normalize_pdb_id, is_tty


@dataclass
class RcsbEndpoints:
    search: str = "https://search.rcsb.org/rcsbsearch/v2/query"
    data_entry: str = "https://data.rcsb.org/rest/v1/core/entry"
    holdings_status: str = "https://data.rcsb.org/rest/v1/holdings/status"
    files_download: str = "https://files.rcsb.org/download"


class RateLimitDecision:
    def __init__(self, action: str, new_max_rps: Optional[float]) -> None:
        self.action = action
        self.new_max_rps = new_max_rps


class Requester:
    def __init__(self, max_rps: float, interactive: bool, warnings: WarningCollector, endpoints: Optional[RcsbEndpoints] = None) -> None:
        self.session = requests.Session()
        self.max_rps = float(max_rps)
        if self.max_rps <= 0:
            self.max_rps = 0.5
        self.interactive = bool(interactive)
        self.warnings = warnings
        if endpoints is None:
            self.endpoints = RcsbEndpoints()
        else:
            self.endpoints = endpoints
        self._last_req_time = 0.0

    def _throttle(self) -> None:
        if self.max_rps <= 0:
            return
        now = time.time()
        min_dt = 1.0 / self.max_rps
        dt = now - self._last_req_time
        if dt < min_dt:
            time.sleep(min_dt - dt)
        self._last_req_time = time.time()

    def _handle_429(self, url: str) -> RateLimitDecision:
        self.warnings.add("RCSB rate limit detected (HTTP 429).", url=url, max_rps=self.max_rps)
        if self.interactive and is_tty():
            print("")
            print("PDBStEx WARNING: Rate limit detected (HTTP 429).")
            print("Options:")
            print("  R  Reduce max requests/second and retry")
            print("  C  Retry without changes")
            print("  S  Skip this request")
            while True:
                ans = input("Choose [R/C/S]: ").strip().upper()
                if ans == "R":
                    new_rps = self.max_rps * 0.5
                    if new_rps < 0.1:
                        new_rps = 0.1
                    return RateLimitDecision(action="reduce", new_max_rps=new_rps)
                if ans == "C":
                    return RateLimitDecision(action="continue", new_max_rps=None)
                if ans == "S":
                    return RateLimitDecision(action="skip", new_max_rps=None)
        new_rps2 = self.max_rps * 0.5
        if new_rps2 < 0.1:
            new_rps2 = 0.1
        return RateLimitDecision(action="reduce", new_max_rps=new_rps2)

    def get_json(self, url: str, timeout: float = 60.0, max_retries: int = 6) -> Optional[dict[str, Any]]:
        attempt = 0
        while attempt <= max_retries:
            self._throttle()
            try:
                resp = self.session.get(url, timeout=timeout, headers={"Accept": "application/json"})
            except Exception as e:
                self.warnings.add("HTTP GET failed.", url=url, error=str(e), attempt=attempt)
                sleep_backoff(attempt)
                attempt += 1
                continue
            if resp.status_code == 200:
                try:
                    return resp.json()
                except Exception as e:
                    self.warnings.add("JSON parsing failed.", url=url, error=str(e))
                    return None
            if resp.status_code == 429:
                decision = self._handle_429(url)
                if decision.action == "skip":
                    return None
                if decision.action == "reduce" and decision.new_max_rps is not None:
                    self.max_rps = float(decision.new_max_rps)
                sleep_backoff(attempt)
                attempt += 1
                continue
            if resp.status_code in [500, 502, 503, 504]:
                self.warnings.add("Server error, retrying.", url=url, status=resp.status_code, attempt=attempt)
                sleep_backoff(attempt)
                attempt += 1
                continue
            if resp.status_code == 404:
                return None
            self.warnings.add("HTTP GET unexpected status.", url=url, status=resp.status_code, body=resp.text[:400])
            return None
        return None

    def post_json(self, url: str, payload: dict[str, Any], timeout: float = 60.0, max_retries: int = 6) -> Optional[dict[str, Any]]:
        attempt = 0
        while attempt <= max_retries:
            self._throttle()
            try:
                resp = self.session.post(url, json=payload, timeout=timeout, headers={"Accept": "application/json"})
            except Exception as e:
                self.warnings.add("HTTP POST failed.", url=url, error=str(e), attempt=attempt)
                sleep_backoff(attempt)
                attempt += 1
                continue
            if resp.status_code == 200:
                try:
                    return resp.json()
                except Exception as e:
                    self.warnings.add("JSON parsing failed.", url=url, error=str(e))
                    return None
            if resp.status_code == 204:
                return {"result_set": [], "total_count": 0}
            if resp.status_code == 429:
                decision = self._handle_429(url)
                if decision.action == "skip":
                    return None
                if decision.action == "reduce" and decision.new_max_rps is not None:
                    self.max_rps = float(decision.new_max_rps)
                sleep_backoff(attempt)
                attempt += 1
                continue
            if resp.status_code in [500, 502, 503, 504]:
                self.warnings.add("Server error, retrying.", url=url, status=resp.status_code, attempt=attempt)
                sleep_backoff(attempt)
                attempt += 1
                continue
            self.warnings.add("HTTP POST unexpected status.", url=url, status=resp.status_code, body=resp.text[:400])
            return None
        return None

    def download_bytes(self, url: str, timeout: float = 120.0, max_retries: int = 6) -> Optional[bytes]:
        attempt = 0
        while attempt <= max_retries:
            self._throttle()
            try:
                resp = self.session.get(url, timeout=timeout, headers={"Accept": "*/*"})
            except Exception as e:
                self.warnings.add("HTTP download failed.", url=url, error=str(e), attempt=attempt)
                sleep_backoff(attempt)
                attempt += 1
                continue
            if resp.status_code == 200:
                return resp.content
            if resp.status_code == 429:
                decision = self._handle_429(url)
                if decision.action == "skip":
                    return None
                if decision.action == "reduce" and decision.new_max_rps is not None:
                    self.max_rps = float(decision.new_max_rps)
                sleep_backoff(attempt)
                attempt += 1
                continue
            if resp.status_code in [500, 502, 503, 504]:
                self.warnings.add("Server error, retrying download.", url=url, status=resp.status_code, attempt=attempt)
                sleep_backoff(attempt)
                attempt += 1
                continue
            if resp.status_code == 404:
                return None
            self.warnings.add("HTTP download unexpected status.", url=url, status=resp.status_code, body=resp.text[:400])
            return None
        return None


def build_search_request(
    terms: list[str],
    methods: list[str],
    max_resolution: float,
    include_no_ligand: bool,
    max_rows: int,
) -> dict[str, Any]:
    nodes: list[dict[str, Any]] = []
    for term in terms:
        t = str(term).strip()
        if t != "":
            nodes.append({"type": "terminal", "service": "full_text", "parameters": {"value": t}})

    if methods is not None and len(methods) > 0:
        m_nodes: list[dict[str, Any]] = []
        for m in methods:
            mv = str(m).strip()
            if mv != "":
                m_nodes.append(
                    {
                        "type": "terminal",
                        "service": "text",
                        "parameters": {
                            "attribute": "rcsb_entry_info.experimental_method",
                            "operator": "exact_match",
                            "value": mv,
                        },
                    }
                )
        if len(m_nodes) == 1:
            nodes.append(m_nodes[0])
        elif len(m_nodes) > 1:
            nodes.append({"type": "group", "logical_operator": "or", "nodes": m_nodes})

    if max_resolution is not None and max_resolution > 0:
        nodes.append(
            {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_entry_info.resolution_combined",
                    "operator": "less_or_equal",
                    "value": float(max_resolution),
                },
            }
        )

    if include_no_ligand is False:
        nodes.append(
            {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_entry_info.nonpolymer_entity_count",
                    "operator": "greater_or_equal",
                    "value": 1,
                },
            }
        )

    query: dict[str, Any]
    if len(nodes) == 0:
        query = {}
    elif len(nodes) == 1:
        query = nodes[0]
    else:
        query = {"type": "group", "logical_operator": "and", "nodes": nodes}

    req: dict[str, Any] = {
        "return_type": "entry",
        "query": query,
        "request_options": {
            "paginate": {"start": 0, "rows": int(max_rows)},
            "results_verbosity": "compact",
        },
    }
    return req


def search_entries(
    requester: Requester,
    terms: list[str],
    methods: list[str],
    max_resolution: float,
    include_no_ligand: bool,
    warnings: WarningCollector,
    max_entries: int,
) -> list[str]:
    url = requester.endpoints.search
    rows = 10000
    if max_entries > 0 and max_entries < rows:
        rows = max_entries

    req = build_search_request(terms=terms, methods=methods, max_resolution=max_resolution, include_no_ligand=include_no_ligand, max_rows=rows)

    resp = requester.post_json(url, req, timeout=90.0)
    if resp is None:
        warnings.add("Search API failed; retrying with a simpler query.")
        req2 = build_search_request(terms=terms, methods=[], max_resolution=0.0, include_no_ligand=True, max_rows=rows)
        resp = requester.post_json(url, req2, timeout=90.0)
        if resp is None:
            warnings.add("Search API failed again; returning empty result set.")
            return []

    out: list[str] = []
    if "result_set" in resp and isinstance(resp["result_set"], list):
        for item in resp["result_set"]:
            if isinstance(item, str):
                pid = normalize_pdb_id(item)
                if pid != "":
                    out.append(pid)
            elif isinstance(item, dict) and "identifier" in item:
                pid2 = normalize_pdb_id(str(item["identifier"]))
                if pid2 != "":
                    out.append(pid2)

    if max_entries > 0:
        out2: list[str] = []
        i = 0
        for x in out:
            out2.append(x)
            i += 1
            if i >= max_entries:
                break
        return out2

    return out


def resolve_obsolete_pdb_id(requester: Requester, pdb_id: str, warnings: WarningCollector) -> str:
    pid = normalize_pdb_id(pdb_id)
    if pid == "":
        return ""
    url = requester.endpoints.holdings_status + "/" + pid
    data = requester.get_json(url, timeout=30.0)
    if data is None:
        return pid

    if isinstance(data, dict):
        status = ""
        if "status" in data:
            status = str(data["status"]).strip().lower()
        if status == "removed" or status == "obsolete":
            replaced = ""
            if "replacedBy" in data and data["replacedBy"] is not None:
                if isinstance(data["replacedBy"], list):
                    if len(data["replacedBy"]) > 0:
                        replaced = str(data["replacedBy"][0]).strip()
                else:
                    replaced = str(data["replacedBy"]).strip()
            if replaced != "":
                rep = normalize_pdb_id(replaced)
                if rep != "":
                    warnings.add("Obsolete PDB ID replaced.", pdb_id=pid, replaced_by=rep)
                    return rep
            warnings.add("PDB ID appears obsolete/removed and no replacement was found.", pdb_id=pid)
    return pid


def fetch_entry_metadata(requester: Requester, pdb_id: str) -> Optional[dict[str, Any]]:
    pid = normalize_pdb_id(pdb_id)
    if pid == "":
        return None
    url = requester.endpoints.data_entry + "/" + pid
    return requester.get_json(url, timeout=90.0)
