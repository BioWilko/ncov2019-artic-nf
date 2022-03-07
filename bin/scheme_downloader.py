#!/usr/bin/env python

import hashlib
import os
import requests
import sys


def check_scheme_hashes(filepath, manifest_hash):
    with open(filepath, "rb") as fh:
        data = fh.read()
        hash_sha256 = hashlib.sha256(data).hexdigest()
    if hash_sha256 != manifest_hash:
        print(
            f"sha256 hash for {str(filepath)} does not match manifest", file=sys.stderr,
        )
        raise SystemExit(1)


def get_scheme(scheme_name, scheme_version="1"):
    """Get and check the ARTIC primer scheme.
    When determining a version, the original behaviour (parsing the scheme_name and
    separating on /V ) is used over a specified scheme_version. If neither are
    provided, the version defaults to 1.
    If 0 is provided as version, the latest scheme will be downloaded.
 
    Parameters
    ----------
    scheme_name : str
        The primer scheme name
    scheme_directory : str
        The directory containing the primer scheme and reference sequence
    scheme_version : str
        The primer scheme version (optional)
    Returns
    -------
    str
        The location of the checked primer scheme
    str
        The location of the checked reference sequence
    str
        The version being used
    """
    # try getting the version from the scheme name (old behaviour)
    if scheme_name.find("/V") != -1:
        scheme_name, scheme_version = scheme_name.split("/V")

    scheme_version = scheme_version.upper().replace("V", "")

    # create the filenames and check they exist
    bed = "scheme/%s/V%s/%s.scheme.bed" % (scheme_name, scheme_version, scheme_name)
    ref = "scheme/%s/V%s/%s.reference.fasta" % (
        scheme_name,
        scheme_version,
        scheme_name,
    )
    if os.path.exists(bed) and os.path.exists(ref):
        return bed, ref, scheme_version

    # if they don't exist, try downloading them to the current directory
    print(
        "could not find primer scheme and reference sequence, downloading",
        file=sys.stderr,
    )

    try:
        manifest = requests.get(
            "https://raw.githubusercontent.com/artic-network/primer-schemes/master/schemes_manifest.json"
        ).json()
    except requests.exceptions.RequestException as error:
        print("Manifest Exception:", error)
        raise SystemExit(2)

    for scheme, scheme_contents in dict(manifest["schemes"]).items():
        if (
            scheme == scheme_name.lower()
            or scheme_name.lower() in scheme_contents["aliases"]
        ):
            print(
                f"\tfound requested scheme:\t{scheme} (using alias {scheme_name})",
                file=sys.stderr,
            )
            if scheme_version == 0:
                print(
                    f"Latest version for scheme {scheme} is -> {scheme_contents['latest_version']}",
                    file=sys.stderr,
                )
                scheme_version = scheme_contents["latest_version"]
            elif scheme_version not in dict(scheme_contents["primer_urls"]).keys():
                print(
                    f"Requested scheme version V{scheme_version} not found; using latest version (V{scheme_contents['latest_version']}) instead",
                    file=sys.stderr,
                )
                scheme_version = scheme_contents["latest_version"]
                bed = "scheme/%s/V%s/%s.scheme.bed" % (
                    scheme_name,
                    scheme_version,
                    scheme_name,
                )
                ref = "scheme/%s/V%s/%s.reference.fasta" % (
                    scheme_name,
                    scheme_version,
                    scheme_name,
                )

            os.makedirs(os.path.dirname(bed), exist_ok=True)
            with requests.get(scheme_contents["primer_urls"][scheme_version]) as fh:
                open(bed, "wt").write(fh.text)

            os.makedirs(os.path.dirname(ref), exist_ok=True)
            with requests.get(scheme_contents["reference_urls"][scheme_version]) as fh:
                open(ref, "wt").write(fh.text)

            check_scheme_hashes(
                bed, scheme_contents["primer_sha256_checksums"][scheme_version]
            )
            check_scheme_hashes(
                ref, scheme_contents["reference_sha256_checksums"][scheme_version]
            )

            return bed, ref, scheme_version

    print(
        f"\tRequested scheme:\t{scheme_name} could not be found, exiting",
        file=sys.stderr,
    )
    raise SystemExit(1)


get_scheme(sys.argv[1], sys.argv[2])
