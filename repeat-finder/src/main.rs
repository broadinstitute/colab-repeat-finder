fn main() {
    for line in std::io::stdin().lines() {
        let line = line.unwrap();
        println!("{line}");
    }
}
