#ifndef INCLUDED_COMPBIO_ABINITIO
#define INCLUDED_COMPBIO_ABINITIO

#include <core/types.hh>
#include <protocols/abinitio/FragmentMover.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pose/Pose.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/WhileMover.hh>
#include <protocols/Protocol.hh>
#include <protocols/evaluation/PoseEvaluator.hh>
#include <protocols/evaluation/PCA.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/silent.fwd.hh>

class compbio_abinitio;
typedef utility::pointer::owning_ptr< compbio_abinitio > compbio_abinitioOP;

class compbio_abinitio : public protocols::Protocol
{
    typedef Protocol Parent;
public:
    static void register_options();
    enum StageID {
        ALL_STAGES = 0,
        STAGE_1,
        STAGE_2,
        STAGE_3a,
        STAGE_3b,
        STAGE_4,
        STAGE_5
    };
    compbio_abinitio (int decoy_num);   
    ~compbio_abinitio();
    void run ();    
 private:
    void setup ();
    void fold (core::pose::Pose &pose);
    bool fold_stage1 (core::pose::Pose &pose);
    bool fold_stage2 (core::pose::Pose &pose);
    bool fold_stage3 (core::pose::Pose &pose);
    bool fold_stage4 (core::pose::Pose &pose);
    bool fold_stage5 (core::pose::Pose &pose);
    void init_on_new_input ();
    void output_models ();
    void generate_extended_pose (core::pose::Pose &pose,
                               std::string const& seq) const;
    void ex_generate_extended_pose (core::pose::Pose &pose,
            core::pose::Pose &oldPose, std::string const& seq) const;
    bool prepare_loop_in_stage3 (core::pose::Pose &pose, Size iteration);
    void recover_low (core::pose::Pose& pose, StageID stage );
    bool contains_stageid (utility::vector1<StageID > vec, StageID query);
    void replace_scorefxn (core::pose::Pose& pose, StageID stage);
    core::scoring::ScoreFunction const& current_scorefxn () const;
    void current_scorefxn (core::scoring::ScoreFunction const&);
    void process_decoy (core::pose::Pose &pose,
                        core::scoring::ScoreFunction const& scorefxn,
                        core::io::silent::SilentStruct &pss) const;
    
    std::string sequence_;                          // target sequence
    core::fragment::FragSetOP fragset_large_;       // probably 9mers
    core::fragment::FragSetOP fragset_small_top25_; // top 25% 3mers
    core::fragment::FragSetOP fragset_small_;       // probably 3mers
    core::fragment::FragSetOP fragset_templates_;   // homolog fragments 
    utility::vector1< StageID > recover_low_stages_;

    protocols::moves::MonteCarloOP mc_;          // Monte Carlo protocol
    core::Real temperature_;                     // temperature

    // No. of cycles to be performed by  Monte Carlo in 5 stages
    core::Size num_cycles_stage1;
    core::Size num_cycles_stage2;
    core::Size num_cycles_stage3;
    core::Size num_cycles_stage4;
    core::Size num_cycles_stage5;

    // Total trials conducted by Monte Carlo
    core::Size total_trials_;
    
    // Scoring functions
    core::scoring::ScoreFunctionOP score_stage1_;  //score0
    core::scoring::ScoreFunctionOP score_stage2_;  //score1
    core::scoring::ScoreFunctionOP score_stage3a_; //score2
    core::scoring::ScoreFunctionOP score_stage3b_; //score5
    core::scoring::ScoreFunctionOP score_stage4_;  //score3
    core::scoring::ScoreFunctionOP score_stage5_;  //score3 

    // Large and small Fragments
    protocols::abinitio::FragmentMoverOP brute_move_small_;
    protocols::abinitio::FragmentMoverOP brute_move_large_;
    protocols::abinitio::FragmentMoverOP smooth_move_small_;

    // Trial movers (to alter/modify protein structure)
    protocols::moves::TrialMoverOP trial_large_;
    protocols::moves::TrialMoverOP trial_small_;
    protocols::moves::TrialMoverOP smooth_trial_small_;

    // switches
    bool apply_large_frags_;
    bool short_insert_region_;
    bool skip_stage1;
    bool skip_stage2;
    bool skip_stage3;
    bool skip_stage4;
    bool skip_stage5;
    bool init_for_input_yet_;
    
    core::scoring::ScoreFunctionOP centroid_scorefunction_; 
    core::scoring::ScoreFunctionOP fullatom_scorefunction_;
    core::kinematics::MoveMapOP movemap;
    core::pose::PoseOP input_pose;            
    core::pose::PoseOP native_pose;
    std::string tag;
    // a bunch of PoseEvaluators for process_decoy() --- if available
	protocols::evaluation::MetaPoseEvaluatorOP evaluator_;
    // a score file ( written to in process_decoy )
	core::io::silent::SilentFileDataOP silent_score_file_;
};

/*----------------------------class hConvergenceCheck--------------------------
 * @brief (helper) functor class which keeps track of old pose for
 * the convergence check in stage3 cycles @detail calls of operator (
 * pose ) compare the class hConvergenceCheck;
 */
class compbio_convergence_check;
typedef  utility::pointer::owning_ptr <compbio_convergence_check>
                                               compbio_convergence_checkOP;
class compbio_convergence_check : public protocols::moves::PoseCondition {
public:
	compbio_convergence_check() : bInit_( false ), ct_( 0 ) {}
	void reset() { ct_ = 0; bInit_ = false; }
	void set_trials(protocols::moves::TrialMoverOP trin ) {
		trials_ = trin;
		runtime_assert(trials_->keep_stats_type() <protocols::moves::no_stats);
		last_move_ = 0;
	}
	virtual bool operator() ( const core::pose::Pose & pose );
private:
	core::pose::Pose very_old_pose_;
	bool bInit_;
	Size ct_;
	protocols::moves::TrialMoverOP trials_;
	Size last_move_;
};

#endif  /*INCLUDED_COMPBIO_ABINITIO*/
