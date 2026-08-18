import pod5

def open_pod5(path: str) -> pod5.Reader:
    return pod5.Reader(path)

def get_signal(reader: pod5.Reader, read_id: str, calibrated: bool = False):
    record = next(
        reader.reads(
            selection=[read_id],
            missing_ok=False,
            preload={"samples"},
        )
    )
    # there was a minknow and pod5 version once, that directly normalized the DACS values instead of the pA
    # for calibrated: check if the shift and scale in the basecalled bams are matching the DACS values or the pA values
    return record.signal_pa if calibrated else record.signal
