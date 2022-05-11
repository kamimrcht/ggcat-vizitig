use crate::execution_manager::executor::Executor;
use crate::execution_manager::executor_address::{ExecutorAddress, WeakExecutorAddress};
use crate::execution_manager::objects_pool::PoolObject;
use crate::execution_manager::packet::{Packet, PacketsPool};
use std::any::TypeId;
use std::borrow::Borrow;
use std::hash::{Hash, Hasher};
use std::ops::Deref;
use std::sync::Arc;

pub type GenericExecutor = Box<dyn ExecutionManagerTrait>;

#[derive(Eq, PartialEq, Debug)]
pub enum ExecutionStatus {
    NoMorePacketsNoProgress,
    NoMorePackets,
    MorePackets,
    OutputPoolFull,
    Finished,
}

pub trait ExecutionManagerTrait: Send + Sync {
    fn execute(&mut self, wait: bool) -> ExecutionStatus;

    fn get_weak_address(&self) -> &WeakExecutorAddress;

    fn get_packets_queue_size(&self) -> usize;

    fn is_finished(&self) -> bool;

    fn exec_type_id(&self) -> TypeId;
}

impl Borrow<WeakExecutorAddress> for dyn ExecutionManagerTrait {
    fn borrow(&self) -> &WeakExecutorAddress {
        self.get_weak_address()
    }
}

impl Hash for dyn ExecutionManagerTrait {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.get_weak_address().hash(state)
    }
}

impl PartialEq for dyn ExecutionManagerTrait {
    fn eq(&self, other: &Self) -> bool {
        self.get_weak_address() == other.get_weak_address()
    }
}

impl Eq for dyn ExecutionManagerTrait {}

impl dyn ExecutionManagerTrait {
    pub fn downcast<E: Executor>(&self) -> &ExecutionManager<E> {
        assert_eq!(self.exec_type_id(), TypeId::of::<E>());
        unsafe { &*(self as *const dyn ExecutionManagerTrait as *const ExecutionManager<E>) }
    }
}

pub struct ExecutionManager<E: Executor> {
    executor: PoolObject<E>,
    weak_address: WeakExecutorAddress,
    pool: Option<Arc<PoolObject<PacketsPool<E::OutputPacket>>>>,
    packets_queue: Box<dyn Fn() -> (Option<Packet<E::InputPacket>>, Option<usize>) + Sync + Send>,
    build_info: Option<E::BuildParams>,
    output_fn: Box<dyn (Fn(ExecutorAddress, Packet<E::OutputPacket>)) + Sync + Send>,
    notify_drop: Box<dyn Fn() + Sync + Send>,
    is_finished: bool,
    queue_size: usize,
}

impl<E: Executor + 'static> ExecutionManager<E> {
    pub fn new(
        executor: PoolObject<E>,
        build_info: E::BuildParams,
        pool: Option<Arc<PoolObject<PacketsPool<E::OutputPacket>>>>,
        packets_queue: Box<
            dyn Fn() -> (Option<Packet<E::InputPacket>>, Option<usize>) + Sync + Send,
        >,
        queue_size: usize,
        weak_address: WeakExecutorAddress,
        output_fn: impl Fn(ExecutorAddress, Packet<E::OutputPacket>) + Sync + Send + 'static,
        notify_drop: impl Fn() + Sync + Send + 'static,
    ) -> GenericExecutor {
        Box::new(Self {
            executor,
            weak_address,
            pool,
            packets_queue,
            build_info: Some(build_info),
            output_fn: Box::new(output_fn),
            notify_drop: Box::new(notify_drop),
            is_finished: false,
            queue_size,
        })
    }
}

impl<E: Executor> ExecutionManagerTrait for ExecutionManager<E> {
    fn execute(&mut self, _wait: bool) -> ExecutionStatus {
        // let available_pool_items = self
        //     .pool
        //     .as_ref()
        //     .map(|p| p.get_available_items())
        //     .unwrap_or(0);
        //
        // while self.executor.required_pool_items() as i64 > available_pool_items {
        //     if wait {
        //         self.pool
        //             .as_ref()
        //             .unwrap()
        //             .wait_for_item_timeout(Duration::from_millis(100));
        //         wait = false;
        //     } else {
        //         return ExecutionStatus::OutputPoolFull;
        //     }
        // }

        if let Some(build_info) = self.build_info.take() {
            self.executor.pre_execute(
                build_info,
                || self.pool.as_ref().unwrap().alloc_packet(false),
                || {
                    self.pool
                        .as_ref()
                        .unwrap()
                        .alloc_packet(E::STRICT_POOL_ALLOC)
                },
                self.output_fn.deref(),
            );
            if self.executor.is_finished() {
                self.is_finished = true;
                return ExecutionStatus::Finished;
            }
        }

        let (packet, queue_size) = (self.packets_queue)();

        let packet_is_none = packet.is_none();

        if let Some(packet) = packet {
            self.executor.execute(
                packet,
                || {
                    self.pool
                        .as_ref()
                        .unwrap()
                        .alloc_packet(E::STRICT_POOL_ALLOC)
                },
                self.output_fn.deref(),
            );
        }

        self.queue_size = queue_size.unwrap_or(0);

        if self.executor.is_finished() || queue_size.is_none() {
            self.is_finished = true;
            ExecutionStatus::Finished
        } else if packet_is_none {
            ExecutionStatus::NoMorePacketsNoProgress
        } else {
            if self.queue_size > 0 {
                ExecutionStatus::MorePackets
            } else {
                ExecutionStatus::NoMorePackets
            }
        }
    }

    fn get_weak_address(&self) -> &WeakExecutorAddress {
        &self.weak_address
    }

    fn is_finished(&self) -> bool {
        self.is_finished || self.executor.is_finished()
    }

    fn get_packets_queue_size(&self) -> usize {
        self.queue_size
    }

    fn exec_type_id(&self) -> TypeId {
        TypeId::of::<E>()
    }
}

impl<E: Executor> Drop for ExecutionManager<E> {
    fn drop(&mut self) {
        self.executor.finalize(self.output_fn.deref());
        (self.notify_drop)();
    }
}
