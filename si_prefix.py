SI_PREFIXES = [
    ("n", 1e-9),
    ("u", 1e-6),
    ("m", 1e-3),
    ("", 1.0),
    ("k", 1e3),
    ("M", 1e6),
    ("G", 1e9),
]


def prefix_labels(unit: str) -> list[str]:
    return [f"{prefix}{unit}" if prefix else unit for prefix, _ in SI_PREFIXES]


def prefix_multiplier(label: str, unit: str) -> float:
    for prefix, multiplier in SI_PREFIXES:
        expected = f"{prefix}{unit}" if prefix else unit
        if label == expected:
            return multiplier
    return 1.0
