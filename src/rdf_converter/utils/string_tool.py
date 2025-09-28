from typing import Optional

def is_not_empty(s: Optional[str]) -> bool:
    return s is not None and s.strip() != ''

