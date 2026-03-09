def set_date_gregorian(year, month, day, hours=0, minutes=0, seconds=0):
    """
    Computes time corresponding to date for Gregorian calendar.
    """

    base_year = 1601

    # Check for valid date and time
    if (
        seconds > 59
        or seconds < 0
        or minutes > 59
        or minutes < 0
        or hours > 23
        or hours < 0
        or day < 1
        or month > 12
        or month < 1
        or year < base_year
    ):

        errstring = f"year,mon,day,hour,min,sec {year} {month} {day} {hours} {minutes} {seconds} not a valid date."
        raise ValueError(errstring)

    # if month != 2 and any([day > month_day for month_day in days_per_month]):
    #    raise ValueError(f"month ({month}) does not have {day} days.")
    if day > days_per_month[month - 1]:
        raise ValueError(f"month ({month}) does not have {day} days.")
    # Check for leap year
    leap = is_leap_year(year)

    if month == 2 and (day > 29 or (not leap and day > 28)):
        raise ValueError(
            f"month ({month}) does not have {day} days in a non-leap year."
        )

    # Compute number of leap years fully past since base_year
    nleapyr = (
        (year - base_year) // 4 - (year - base_year) // 100 + (year - base_year) // 400
    )

    # Count up days in this year
    ndays = sum(
        days_per_month[m - 1] + (1 if leap and m == 2 else 0) for m in range(1, month)
    )

    totseconds = seconds + 60 * (minutes + 60 * (hours))
    totdays = day - 1 + ndays + 365 * (year - base_year - nleapyr) + 366 * nleapyr

    return totseconds, totdays


def is_leap_year(year):
    """
    Checks if the given year is a leap year.
    """
    if year % 4 != 0:
        return False
    elif year % 100 != 0:
        return True
    elif year % 400 != 0:
        return False
    else:
        return True


days_per_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]  # Days in each month
