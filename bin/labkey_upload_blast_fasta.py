#!/usr/bin/env python3
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "labkey>=4.0.1",
# ]
# ///

import argparse
import csv
import os
import sys
from datetime import datetime

from py_nvd.labkey_io import PartialInsertError, insert_records, rows_present


def combo_already_uploaded(
    query_api,
    schema,
    table,
    experiment,
    sample_id,
    query_class,
) -> bool:
    """True if the FASTA list already holds this unit's rows (its ledger).

    The unit is (experiment, sample_id, query_class), matching the BLAST hits
    upload. Keying on (experiment, sample_id) alone was correct only while the
    list held contigs: now that every queried class is uploaded, the contig
    batch would make each later class for that sample look already uploaded and
    its rows would be dropped with no error. Reuses the guard script's filter
    fallbacks (QueryFilter objects, then filter-string syntax) so this still
    works against older labkey-api-python versions.
    """
    return rows_present(
        query_api,
        schema,
        table,
        {
            "experiment": experiment,
            "sample_id": sample_id,
            "query_class": query_class,
        },
    )


def _write_log(log_entries: list[str]) -> None:
    """Write the accumulated log entries to the upload log file and stdout."""
    text = "\n".join(log_entries)
    with open("fasta_labkey_upload.log", "w") as f:
        f.write(text)
    print(text)


def main():
    parser = argparse.ArgumentParser(description="Upload FASTA CSVs to LabKey.")
    parser.add_argument("--experiment-id", required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--query-class", required=True)
    parser.add_argument("--labkey-server", required=True)
    parser.add_argument("--labkey-project-name", required=True)
    parser.add_argument("--labkey-api-key", required=True)
    parser.add_argument("--labkey-schema", required=True)
    parser.add_argument("--table-name", default="fasta_hits_test")
    parser.add_argument(
        "--insert-batch-size",
        type=int,
        default=1000,
        help=(
            "Rows per LabKey insert call. Read-derived query classes produce "
            "payloads large enough that a single call can hang the server."
        ),
    )
    args = parser.parse_args()

    log_entries = [
        f"LabKey FASTA Upload Log - {datetime.now()}",
        f"Experiment ID: {args.experiment_id}",
        f"Sample: {args.sample_id}",
        f"Query class: {args.query_class}",
        f"Server: {args.labkey_server}",
        f"Project: {args.labkey_project_name}",
        f"Target Table: {args.table_name}",
        "=" * 80,
    ]

    upload_enabled = bool(
        args.labkey_server and args.labkey_project_name and args.labkey_api_key,
    )
    if upload_enabled:
        try:
            from labkey.api_wrapper import APIWrapper

            log_entries.append("LabKey API wrapper found - attempting real upload")
            api = APIWrapper(
                args.labkey_server,
                args.labkey_project_name,
                api_key=args.labkey_api_key,
            )
        except ImportError as e:
            log_entries.append(f"ERROR: LabKey API wrapper not available - {e!s}")
            log_entries.append("Falling back to simulation mode")
            upload_enabled = False
    else:
        log_entries.append(
            "No LabKey credentials provided - running in simulation mode",
        )

    # The destination list is its own completion ledger. If this
    # (sample_id, query_class) unit already has rows there, a prior run
    # uploaded it: skip re-inserting rather than duplicate the FASTA list.
    if upload_enabled and combo_already_uploaded(
        api.query,
        args.labkey_schema,
        args.table_name,
        int(args.experiment_id),
        args.sample_id,
        args.query_class,
    ):
        log_entries.append(
            f"SKIP: combo already uploaded (exp={args.experiment_id}, "
            f"sample={args.sample_id}, query_class={args.query_class}); no insert.",
        )
        _write_log(log_entries)
        return

    csv_files = [
        f for f in os.listdir(".") if f.endswith(".csv") and "fasta" in f.lower()
    ]
    log_entries.append(f"Found {len(csv_files)} FASTA CSV files: {csv_files}")

    total_records_processed = 0
    total_records_uploaded = 0

    for csv_file in sorted(csv_files):
        log_entries.append(f"\nProcessing FASTA file: {csv_file}")
        record_count = 0

        if os.path.getsize(csv_file) > 0:
            with open(csv_file) as f:
                reader = csv.DictReader(f)
                records = list(reader)

                record_count = len(records)
                total_records_processed += record_count

                if record_count > 0:
                    log_entries.append(f"  Records: {record_count}")
                    log_entries.append(
                        f"  Sample fields: {', '.join(list(records[0].keys())[:5])}",
                    )

                    if upload_enabled:
                        # Sent in batches: read-derived query classes make
                        # these payloads large enough that a single insert call
                        # can hang the server. Batching gives up atomicity, so a
                        # mid-way failure is reported loudly rather than rolled
                        # back.
                        try:
                            insert_records(
                                api.query,
                                args.labkey_schema,
                                args.table_name,
                                records,
                                batch_size=args.insert_batch_size,
                            )
                        except PartialInsertError as e:
                            log_entries.append(f"  Upload: ERROR - {e!s}")
                            log_entries.append(
                                f"\nFASTA UPLOAD FAILED - {e.rows_committed} of "
                                f"{e.total_rows} rows were COMMITTED and were NOT "
                                f"rolled back.\n"
                                f"The destination list is keyed on row presence, "
                                f"so a retry will treat "
                                f"experiment={args.experiment_id} "
                                f"sample={args.sample_id} "
                                f"query_class={args.query_class} as already "
                                f"uploaded and "
                                f"silently skip the remaining "
                                f"{e.total_rows - e.rows_committed} rows.\n"
                                f"Delete that unit's rows from "
                                f"{args.labkey_schema}.{args.table_name} before "
                                f"re-running.",
                            )
                            _write_log(log_entries)
                            print(
                                f"ERROR: LabKey insert failed for "
                                f"sample={args.sample_id} query_class={args.query_class} "
                                f"after committing "
                                f"{e.rows_committed}/{e.total_rows} rows; delete "
                                f"this unit's rows before retrying: {e!s}",
                                file=sys.stderr,
                            )
                            sys.exit(1)
                        except Exception as e:
                            log_entries.append(f"  Upload: ERROR - {e!s}")
                            log_entries.append(
                                "\nFASTA UPLOAD FAILED - insert error, see above (no rows committed)",
                            )
                            _write_log(log_entries)
                            print(
                                f"ERROR: LabKey insert failed for "
                                f"sample={args.sample_id}: {e!s}",
                                file=sys.stderr,
                            )
                            sys.exit(1)
                        log_entries.append(
                            f"  Upload: SUCCESS ({len(records)} records)"
                        )
                        total_records_uploaded += record_count

                    else:
                        num_batches = (record_count + 999) // 1000
                        log_entries.append(
                            f"  Would upload in {num_batches} batch(es) (SIMULATION)",
                        )
                        total_records_uploaded += record_count
                else:
                    log_entries.append("  No records found in file")
        else:
            log_entries.append("  Empty file - no records to upload")

    log_entries += [
        "\n" + "=" * 80,
        "FASTA DATA UPLOAD SUMMARY",
        f"Files processed: {len(csv_files)}",
        f"Total records processed: {total_records_processed}",
        f"Total records uploaded: {total_records_uploaded}"
        if upload_enabled
        else f"Total records that would be uploaded: {total_records_uploaded}",
        "FASTA UPLOAD COMPLETE"
        if upload_enabled
        else "FASTA SIMULATION COMPLETE - No actual upload performed",
    ]

    _write_log(log_entries)


if __name__ == "__main__":
    main()
