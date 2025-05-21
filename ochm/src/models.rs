pub(crate) mod conv_lstm {
    use burn::config::Config;
    use burn::module::Module;
    use burn::nn::conv::{Conv2d, Conv2dConfig};
    use burn::nn::{
        BatchNorm, BatchNormConfig, Linear, LinearConfig, Lstm, LstmConfig,
        Relu,
    };
    use burn::prelude::Backend;
    use burn::tensor::Tensor;

    #[derive(Module, Debug)]
    pub(crate) struct ConvLstmModel<B: Backend> {
        relu: Relu,
        conv1: Conv2d<B>,
        bn1: BatchNorm<B, 2>,
        conv2: Conv2d<B>,
        bn2: BatchNorm<B, 2>,
        conv3: Conv2d<B>,
        bn3: BatchNorm<B, 2>,
        lstm: Lstm<B>,
        linear: Linear<B>,
        class_weights: Option<Vec<f32>>,
    }

    impl<B: Backend> ConvLstmModel<B> {
        pub(crate) fn forward(&self, pileups: Tensor<B, 4>) -> Tensor<B, 2> {
            let x = self.conv1.forward(pileups);
            let x = self.bn1.forward(x);
            let x = self.relu.forward(x);

            let x = self.conv2.forward(x);
            let x = self.bn2.forward(x);
            let x = self.relu.forward(x);

            let x = self.conv3.forward(x);
            let x = self.bn3.forward(x);
            let x = self.relu.forward(x);

            let x = x.permute([0, 3, 2, 1]).squeeze::<3>(1);
            let (x, _) = self.lstm.forward(x, None);
            let [batch_size, sequence_length, hidden_size] = x.dims();
            let x = x
                .slice([0..batch_size, (sequence_length - 1)..sequence_length])
                .squeeze::<2>(1);
            let x = x.reshape([batch_size, hidden_size]);
            let x = self.linear.forward(x);
            x
        }
    }

    #[derive(Config, Debug)]
    pub struct ConvLstmModelConfig {
        num_features: usize,
        hidden_size: usize,
        num_classes: usize,
    }

    impl ConvLstmModelConfig {
        pub(crate) fn init<B: Backend>(
            &self,
            class_weights: Option<Vec<f32>>,
            device: &B::Device,
        ) -> ConvLstmModel<B> {
            let conv1 = Conv2dConfig::new([1, 32], [1, self.num_features])
                .with_stride([1, 1])
                .init::<B>(&device);
            let bn1 = BatchNormConfig::new(32).init::<B, 2>(&device);
            let conv2 = Conv2dConfig::new([32, 128], [4, 1])
                .with_stride([1, 1])
                .init::<B>(&device);
            let bn2 = BatchNormConfig::new(128).init::<B, 2>(&device);

            let conv3 = Conv2dConfig::new([128, 128], [4, 1])
                .with_stride([1, 1])
                .init::<B>(&device);
            let bn3 = BatchNormConfig::new(128).init::<B, 2>(&device);

            let lstm = LstmConfig::new(128, self.hidden_size, false)
                .init::<B>(&device);
            let linear =
                LinearConfig::new(self.hidden_size, 2).init::<B>(&device);

            ConvLstmModel {
                relu: Relu::new(),
                conv1,
                conv2,
                conv3,
                bn1,
                bn2,
                bn3,
                lstm,
                linear,
                class_weights,
            }
        }
    }

    #[cfg(test)]
    #[cfg(feature = "tch")]
    mod conv_lstm_tests {
        use burn::backend::libtorch::LibTorchDevice;
        use burn::backend::LibTorch;
        use burn::module::Module;
        use burn::record::{CompactRecorder, Recorder};

        use crate::models::conv_lstm::{ConvLstmModel, ConvLstmModelConfig};
        use crate::util::ModelConfiguration;

        #[test]
        fn test_load_conv_lstm_model() {
            let device = LibTorchDevice::Cpu;
            let model_weights = CompactRecorder::new()
                .load(
                    std::path::Path::new("/tmp/scratch/model_updated.mpk")
                        .to_path_buf(),
                    &device,
                )
                .unwrap();
            let model_config =
                ModelConfiguration::from_path("/tmp/scratch/model_config.json")
                    .unwrap();
            let model: ConvLstmModel<LibTorch> = ConvLstmModelConfig::new(
                model_config.num_features,
                model_config.hidden_size,
                model_config.num_classes,
            )
            .init(None, &device)
            .load_record(model_weights);
        }
    }
}
