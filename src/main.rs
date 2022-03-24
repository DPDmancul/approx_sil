/*
 *  Copyright 2022 Davide Peressoni
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *    http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */

use rand_distr::{Bernoulli, Distribution};

pub struct Clustering<const N: usize>(pub Vec<Cluster<N>>);
pub struct Cluster<const N: usize>(Vec<Point<N>>, Vec<(Point<N>, f64)>);

#[derive(Copy, Clone)]
pub struct Point<const N: usize>(pub [f64; N]);

impl<const N: usize> Clustering<N> {
    pub fn new<T: Into<Point<N>>>(v: Vec<Vec<T>>, t: u32, δ: f64) -> Self {
        let k = v.len() as u32;
        Self(v.into_iter().map(|c| Cluster::new(c, k, t, δ)).collect())
    }

    pub fn sil(&self) -> f64 {
        let mut i = 0;
        self.0
            .iter()
            .map(|c| {
                i += 1;
                c.sum_sil(self.0[..i - 1].iter().chain(self.0[i..].iter()).collect())
            })
            .sum::<f64>()
            / self.0.iter().map(|c| c.0.len()).sum::<usize>() as f64
    }
}

impl<const N: usize> Cluster<N> {
    pub fn new<T: Into<Point<N>>>(c: Vec<T>, k: u32, t: u32, δ: f64) -> Self {
        let c: Vec<_> = c.into_iter().map(Into::into).collect();
        let sample = if c.len() as u32 <= t {
            c.iter().copied().map(|p| (p, 1.)).collect()
        } else {
            let unif = 1. / c.len() as f64;
            let prob = Bernoulli::new(2. * unif * (2. * k as f64 / δ as f64).ln()).unwrap();
            let s_0: Vec<_> = c
                .iter()
                .filter(|_| prob.sample(&mut rand::thread_rng()))
                .collect();
            let w_0: Vec<_> = s_0
                .iter()
                .map(|p| c.iter().map(|e| p.d(e)).sum::<f64>())
                .collect();
            let prob = c.iter().map(|p| {
                1f64.min(
                    t as f64
                        * unif.max(
                            //γ
                            s_0.iter()
                                .zip(w_0.iter())
                                .map(|(e, w)| p.d(e) / w)
                                .reduce(f64::max)
                                .unwrap(),
                        ),
                )
            });
            c.iter()
                .copied()
                .zip(prob)
                .filter(|(_, p)| Bernoulli::new(*p).unwrap().sample(&mut rand::thread_rng()))
                .collect()
        };
        Self(c, sample)
    }

    fn sum_sil(&self, others: Vec<&Cluster<N>>) -> f64 {
        self.0
            .iter()
            .map(|e| e.sil(self, others.iter().copied()))
            .sum()
    }

    fn w(&self, p: &Point<N>) -> f64 {
        self.1.iter().map(|(e, p_e)| p.d(e) / p_e).sum()
    }
}

impl<const N: usize> From<[f64; N]> for Point<N> {
    fn from(a: [f64; N]) -> Self {
        Self(a)
    }
}

impl<const N: usize> Point<N> {
    fn sil<'a, I: Iterator<Item = &'a Cluster<N>>>(&self, cluster: &Cluster<N>, others: I) -> f64 {
        if cluster.0.len() <= 1 {
            0.
        } else {
            let a = cluster.w(self) / (cluster.0.len() - 1) as f64;
            let b = others.fold(f64::MAX, |min, c| min.min(c.w(self) / c.0.len() as f64));
            (b - a) / a.max(b)
        }
    }

    pub fn d(&self, p: &Point<N>) -> f64 {
        self.0
            .iter()
            .zip(p.0.iter())
            .map(|(a, b)| a - b)
            .map(|x| x * x)
            .sum::<f64>()
            .sqrt()
    }
}

fn main() {
    let c = Clustering::new(
        // exact silhouette 0.19340581656247494
        vec![
            vec![
                [1., 1.],
                [2., 2.],
                [3., 2.],
                [4., 2.],
                [-1., 3.],
                [5., -1.],
                [60.3, 24.5],
                [1., 0.5],
                [5., 3.],
                [20., 21.],
            ],
            vec![
                [-1., -1.],
                [-2., -2.],
                [-3., -2.],
                [-4., -2.],
                [1., -3.],
                [-5., 1.],
                [-60.3, -24.5],
                [-1., -0.5],
                [-5., -3.],
                [-20., -21.],
            ],
        ],
        4, //1024,
        0.1,
    );

    println!("{}", c.sil());
}
