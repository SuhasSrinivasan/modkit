use crate::features::CountsFeatures;
use burn::data::dataloader::batcher::Batcher;
use burn::prelude::Backend;
use burn::tensor::{Shape, Tensor, TensorData};
use derive_new::new;

#[derive(new, Debug, Clone)]
pub(crate) struct PileupBatcher<B: Backend> {
    device: B::Device,
}

#[derive(Debug, Clone)]
pub(crate) struct InfCountsBatch<B: Backend> {
    pub samples: Tensor<B, 3>,
}

impl<B: Backend> Batcher<B, CountsFeatures<()>, InfCountsBatch<B>>
    for PileupBatcher<B>
{
    fn batch(
        &self,
        items: Vec<CountsFeatures<()>>,
        _device: &B::Device,
    ) -> InfCountsBatch<B> {
        let samples =
            items
                .into_iter()
                .map(|xs| {
                    let sh = Shape::new([xs.width, xs.n_features]);
                    let d = TensorData::new(xs.counts, sh);
                    let t = Tensor::<B, 2>::from_data(d, &self.device)
                        .reshape([1, xs.width, xs.n_features]);
                    t
                })
                .collect::<Vec<Tensor<B, 3>>>();
        let samples = Tensor::cat(samples, 0).to_device(&self.device);

        InfCountsBatch { samples }
    }
}
