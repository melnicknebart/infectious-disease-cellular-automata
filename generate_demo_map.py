from pathlib import Path


def main():
    output = Path("data/demo_map.csv")
    output.parent.mkdir(parents=True, exist_ok=True)

    with output.open("w", encoding="utf-8") as f:
        for x in range(150):
            for y in range(150):
                f.write(f"{x},{y}\n")

    print(f"Created synthetic demo map at {output}")


if __name__ == "__main__":
    main()
