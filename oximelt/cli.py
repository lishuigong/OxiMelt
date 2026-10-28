from __future__ import annotations
import argparse
import json
from .core import predict_formula, predict_csv

def main():
    p = argparse.ArgumentParser(description="OxiMelt — oxide melting-point calculator")
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("predict", help="Predict one oxide formula")
    s.add_argument("formula")

    b = sub.add_parser("batch", help="Predict a CSV containing a 'formula' column")
    b.add_argument("input_csv")
    b.add_argument("output_csv")

    args = p.parse_args()
    if args.cmd == "predict":
        r = predict_formula(args.formula)
        print(json.dumps(r.to_flat_dict(), indent=2, ensure_ascii=False))
    else:
        n = predict_csv(args.input_csv, args.output_csv)
        print(f"Processed {n} formula(s) -> {args.output_csv}")

if __name__ == "__main__":
    main()
