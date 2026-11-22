//! Number formatting for display.

/// Fixed decimals with thousands separators for the integer part.
pub fn number(value: f64, decimals: usize) -> String {
    let s = format!("{value:.decimals$}");
    let (int, frac) = s.split_once('.').map_or((s.as_str(), None), |(i, f)| (i, Some(f)));
    let (sign, digits) = int.strip_prefix('-').map_or(("", int), |d| ("-", d));
    let mut grouped = String::new();
    for (i, c) in digits.chars().enumerate() {
        if i > 0 && (digits.len() - i) % 3 == 0 {
            grouped.push(',');
        }
        grouped.push(c);
    }
    match frac {
        Some(f) => format!("{sign}{grouped}.{f}"),
        None => format!("{sign}{grouped}"),
    }
}

/// Percentage with one decimal.
pub fn percent(fraction: f64) -> String {
    format!("{:.1}", fraction * 100.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn groups_thousands() {
        assert_eq!(number(8564.84, 1), "8,564.8");
        assert_eq!(number(-1_234_567.0, 0), "-1,234,567");
        assert_eq!(number(12.0, 2), "12.00");
    }
}
