class SmoothingError(Exception):
    def __init__(self, *args, **kwargs):
        default_message = "The smoothing degree is too high. Retry with a smaller one."

        if args or kwargs:
            super().__init__(*args, **kwargs)
        else:
            super().__init__(default_message)
