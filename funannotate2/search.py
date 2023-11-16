import pyhmmer


def pyhmmer_version():
    return pyhmmer.__version__


def hmmer_search(hmmfile, sequences):
    hmm = next(pyhmmer.plan7.HMMFile(hmmfile))
    results = []
    for top_hits in pyhmmer.hmmsearch([hmm], sequences):
        for hit in top_hits:
            cog = hit.best_domain.alignment.hmm_name.decode()
            domains = []
            for h in hit.domains:
                domains.append(
                    {
                        "hmm_from": h.alignment.hmm_from,
                        "hmm_to": h.alignment.hmm_to,
                        "env_from": h.env_from,
                        "env_to": h.env_to,
                        "score": h.score,
                    }
                )
            results.append(
                {
                    "name": cog,
                    "hit": hit.name.decode(),
                    "bitscore": hit.score,
                    "evalue": hit.evalue,
                    "domains": domains,
                }
            )
    return results
