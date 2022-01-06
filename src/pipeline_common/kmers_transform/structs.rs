use crate::assemble_pipeline::parallel_kmers_merge::KmersFlags;
use crate::colors::colors_manager::ColorsManager;
use crate::config::{BucketIndexType, SortingHashType};
use crate::io::concurrent::intermediate_storage::{IntermediateReadsReader, SequenceExtraData};
use crate::io::concurrent::intermediate_storage_single::IntermediateSequencesStorageSingleBucket;
use crate::io::varint::decode_varint_flags;
use crate::io::{DataReader, DataWriter};
use crate::pipeline_common::kmers_transform::MERGE_BUCKETS_COUNT;
use crate::{CompressedRead, KEEP_FILES};
use crossbeam::queue::SegQueue;
use parallel_processor::lock_free_binary_writer::LockFreeBinaryWriter;
use parallel_processor::memory_fs::file::internal::MemoryFileMode;
use parallel_processor::memory_fs::file::writer::FileWriter;
use parallel_processor::multi_thread_buckets::MultiThreadBuckets;
use parking_lot::Mutex;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Arc;

#[derive(Copy, Clone, Ord, PartialOrd, Eq, PartialEq, Debug)]
pub struct ReadRef {
    pub read_start: usize,
    pub hash: SortingHashType,
}

impl ReadRef {
    #[inline(always)]
    pub fn unpack<'a, 'b, E: SequenceExtraData>(
        &'a self,
        memory: &'b[u8],
        flags_count: usize,
    ) -> (u8, CompressedRead<'b>, E) {
        let mut read_start = self.read_start;
        let (read_len, flags) = decode_varint_flags::<_>(
            || {
                let x = unsafe { memory[read_start] };
                read_start += 1;
                Some(x)
            },
            flags_count,
        )
        .unwrap();

        let read_bases_start = read_start;
        let read_len = read_len as usize;
        let read_len_bytes = (read_len + 3) / 4;

        let read_slice = &memory[read_bases_start..(read_bases_start + read_len_bytes)];

        let read = CompressedRead::new_from_compressed(read_slice, read_len);

        let extra_data = unsafe {
            let extra_slice_ptr = memory.as_ptr().add(read_bases_start + read_len_bytes);
            E::decode_from_pointer(extra_slice_ptr).unwrap()
        };

        (flags, read, extra_data)
    }
}

unsafe impl Sync for ReadRef {}
unsafe impl Send for ReadRef {}

pub struct BucketProcessData<E: SequenceExtraData, R: DataReader> {
    pub reader: IntermediateReadsReader<E, R>,
    pub buckets: MultiThreadBuckets<LockFreeBinaryWriter>,
    vecs_queue: Arc<SegQueue<PathBuf>>,
}

static QUEUE_IDENTIFIER: AtomicU64 = AtomicU64::new(0);

impl<E: SequenceExtraData, R: DataReader> BucketProcessData<E, R> {
    pub fn new(path: impl AsRef<Path>, vecs_queue: Arc<SegQueue<PathBuf>>) -> Self {
        let tmp_dir = path.as_ref().parent().unwrap_or(Path::new("."));
        Self {
            reader: IntermediateReadsReader::<E, R>::new(
                &path,
                !KEEP_FILES.load(Ordering::Relaxed),
            ),
            buckets: MultiThreadBuckets::new(
                MERGE_BUCKETS_COUNT,
                &(
                    PathBuf::from(tmp_dir).join(format!(
                        "vec{}",
                        QUEUE_IDENTIFIER.fetch_add(1, Ordering::Relaxed)
                    )),
                    MemoryFileMode::AlwaysMemory,
                ),
                None,
            ),
            vecs_queue,
        }
    }
}

impl<E: SequenceExtraData, R: DataReader> Drop for BucketProcessData<E, R> {
    fn drop(&mut self) {
        for path in self.buckets.finalize() {
            self.vecs_queue.push(path);
        }
    }
}
