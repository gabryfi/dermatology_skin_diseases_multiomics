process validateReferences {
    tag "reference_validation"

    input:
    path genome
    path annotation
    path checksum_file
    output:
    path "reference_validation.txt"

    script:
    """
    set -euo pipefail

    echo "Reference validation" > reference_validation.txt
    echo "====================" >> reference_validation.txt
    echo "" >> reference_validation.txt

    echo "Genome:" >> reference_validation.txt
    sha256sum ${genome} >> reference_validation.txt

    echo "Annotation:" >> reference_validation.txt
    sha256sum ${annotation} >> reference_validation.txt

    echo "" >> reference_validation.txt
    echo "Checking SHA-256 checksums..." >> reference_validation.txt

    sha256sum -c ${checksum_file}

    echo "" >> reference_validation.txt
    echo "STATUS: PASS" >> reference_validation.txt
    """
}