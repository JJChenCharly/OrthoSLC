# AGENTS.md: `src/TK_SNPmat.py` by L2 norm calculation.

# Current method
- Current all-pairwise SNP count by `src/TK_SNPmat.py` is based on `from itertools import combinations` and a for loop as underlying basic logic.

# Modification to do and how
- NOTE! For following calculation, we will adopt numpy for calculation as package torch is too heavy and not convinient to install. Also keep the multiprocessing strategy to do such calculation for multiple files at same time.
- for one file already aligned by mafft or kalign:
  * I want to transform my input sequences into one-hot like encoded manner like following:
    ```python
    conv_d = {'A': [1., 0., 0., 0.],
              'T': [0., 1., 0., 0.],
              'C': [0., 0., 1., 0.],
              'G': [0., 0., 0., 1.],
              'Y': [0., 0.5, 0.5, 0.], # T, C
              'R': [0.5, 0., 0., 0.5], # A, G
              'S': [0., 0., 0.5, 0.5], # C, G
              'W': [0.5, 0.5, 0., 0.], # A, T
              'K': [0., 0.5, 0., 0.5], # G, T
              'M': [0.5, 0., 0.5, 0.], # A, C
              'B': [0., 1/3, 1/3, 1/3], # T, C, G
              'D': [1/3, 1/3, 0., 1/3], # A, T, G
              'H': [1/3, 1/3, 1/3, 0.], # A, T, C
              'V': [1/3, 0., 1/3, 1/3], # A, C, G
              'N': [0.25, 0.25, 0.25, 0.25],
              '-': [0., 0., 0., 0.]
            }

    pt_path = 'pt/'
    aligned_path = 'aligned/'

    def one_hot_encoding(sequence, conv_d_):
        
        seq_tensor = [conv_d_[base] for base in sequence]
        seq_tensor = torch.tensor(seq_tensor).t()
        
        return seq_tensor

    def DNA_to_pt(in_lst, convert_dict):
        
        for f in in_lst:
            
            the_fasta = SeqIO.to_dict(SeqIO.parse(aligned_path + f, 'fasta'))

            for x in cluster_lst:
              s_seq = str(the_fasta[x].seq)
              b_seq = one_hot_encoding(s_seq, convert_dict)
    ```
  * I will later flatten the tensor of each sequence. For example:
    ```python
    _3_bp_tensor = torch.tensor([
        [1, 2, 3],
        [4, 5, 6],
        [7, 8, 9],
        [10, 11, 12]
    ])

    flattened_3_bp_tensor = _3_bp_tensor.view(1, -1)
    ```
  * Then I will calculate all pairwise distance with something like below:
    ```python
    class VectorDataset(Dataset):
        def __init__(self, vectors, batch_size=32):
            """
            Args:
                vectors: Tensor of shape (total_vectors, vector_dim)
                batch_size: 每个batch包含的向量数量
            """
            self.vectors = vectors
            self.batch_size = batch_size
            self.num_batches = (len(vectors) + batch_size - 1) // batch_size
        
        def __len__(self):
            return self.num_batches
        
        def __getitem__(self, idx):
            start_idx = idx * self.batch_size
            end_idx = min((idx + 1) * self.batch_size, len(self.vectors))
            return self.vectors[start_idx:end_idx]

    class BatchPairwiseDistance(nn.Module):
        def __init__(self):
            super().__init__()
        
        def forward(self, x):
            """
            计算batch内向量的两两距离
            
            Args:
                x: (batch_size, num_vectors, vector_dim)
            """
            return torch.cdist(x, x, p=2)
    ```
    * This will allow faster calculation.