//! Disjoint-set forest over string keys.

use std::collections::{BTreeMap, HashMap};

#[derive(Debug, Default)]
pub struct UnionFind {
    index: HashMap<String, usize>,
    names: Vec<String>,
    parent: Vec<usize>,
    rank: Vec<u8>,
}

impl UnionFind {
    pub fn new() -> Self {
        Self::default()
    }

    fn id(&mut self, key: &str) -> usize {
        if let Some(&i) = self.index.get(key) {
            return i;
        }
        let i = self.parent.len();
        self.index.insert(key.to_owned(), i);
        self.names.push(key.to_owned());
        self.parent.push(i);
        self.rank.push(0);
        i
    }

    fn find_idx(&mut self, mut i: usize) -> usize {
        while self.parent[i] != i {
            self.parent[i] = self.parent[self.parent[i]];
            i = self.parent[i];
        }
        i
    }

    /// Ensures `key` exists as a singleton if it has not been seen.
    pub fn add(&mut self, key: &str) {
        self.id(key);
    }

    pub fn union(&mut self, a: &str, b: &str) {
        let (ia, ib) = (self.id(a), self.id(b));
        let (ra, rb) = (self.find_idx(ia), self.find_idx(ib));
        if ra == rb {
            return;
        }
        match self.rank[ra].cmp(&self.rank[rb]) {
            std::cmp::Ordering::Less => self.parent[ra] = rb,
            std::cmp::Ordering::Greater => self.parent[rb] = ra,
            std::cmp::Ordering::Equal => {
                self.parent[rb] = ra;
                self.rank[ra] += 1;
            }
        }
    }

    /// Maps every key to a canonical component name: the lexicographically
    /// smallest key in its component. Independent of insertion order.
    pub fn components(&mut self) -> BTreeMap<String, String> {
        let n = self.parent.len();
        let mut smallest: HashMap<usize, String> = HashMap::new();
        for i in 0..n {
            let r = self.find_idx(i);
            let name = &self.names[i];
            smallest
                .entry(r)
                .and_modify(|s| {
                    if name < s {
                        *s = name.clone();
                    }
                })
                .or_insert_with(|| name.clone());
        }
        let mut out = BTreeMap::new();
        for i in 0..n {
            let r = self.find_idx(i);
            out.insert(self.names[i].clone(), smallest[&r].clone());
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn merges_transitively() {
        let mut uf = UnionFind::new();
        uf.union("c", "b");
        uf.union("b", "a");
        uf.add("z");
        let comp = uf.components();
        assert_eq!(comp["a"], "a");
        assert_eq!(comp["b"], "a");
        assert_eq!(comp["c"], "a");
        assert_eq!(comp["z"], "z");
    }
}
