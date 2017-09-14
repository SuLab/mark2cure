var channel = Backbone.Radio.channel('mark2cure');

//
// Views
//


RETrainingActionLogic = Object({

  //-- Modules
  0: {
    1:{
      //-- Steps
      onStartUpLogic: function() {
        $('.holder').remove();
        var pos_1 = [$(".paragraph-content").offset(), $(".paragraph-content").height()/2, $(".paragraph-content").width()/2];
        var pos_2 = [$("#c2 .text").offset(), -12, -6];
        var pos_3 = [$("#c1 .concept").offset(), 10, 10];
        var pos_4 = [$("ul a.list-group-item:nth(1)").offset(), -4, -4];
        var pos_5 = [$("#submit_button").offset(), 10, 30];

        _.each([pos_1, pos_2, pos_3, pos_4, pos_5], function(pos, idx) {
          $("body").append("<div class='holder' data-idx='"+idx+"' style='top:"+(pos[0].top+pos[1])+"px; left:"+(pos[0].left+pos[2])+"px;'><p>"+(idx+1)+"</p><div class='dot'></div><div class='pulse'></div></div>");
        });
        $('.holder').on('click mouseenter', function() {
          channel.trigger('training:set:instruction', +$(this).data('idx'));
        });
      },

      onDestroy: function() {
        $('.holder').remove();
      }
    }
  },
  1: {
    //-- Steps
    0: {
      //-- Training
      1: {

        incorrectClick: function(reconcept) {
          this.better_answer();
          this.event_counter++;
        },

        clearClick: function() {
          if(this.event_counter > 0) {
            $('.popover').popover('dispose');
            $relates_el = $('ul a.list-group-item:nth(0)');
            $relates_el.popover({
              container: 'body',
              html: true,
              animation: false,
              content: function() {
                return 'Once again, click on the \'relates to\' box';
              },
              placement: 'left',
            });
            $relates_el.popover('show');
          }
          this.event_counter++;
        },

        choiceClick: function(rechoice_model) {
          if(rechoice_model.get('id') == '8qota4u8hwtcyp65kz9zm0vjyuxwjt12sko084sn') {
            //-- The first "relates to"
            $('.popover').popover('dispose');
            $treats_el = $('ul a.list-group-item:nth(1)');
            $treats_el.popover({
              container: 'body',
              html: true,
              animation: false,
              content: function() {
                return 'When you make a selection, more options may become available. Select the most detailed relationship you can without guessing. <br /><strong>Select \'(may) treat(s)\' to continue</strong>';
              },
              placement: 'bottom'
            });
            $treats_el.popover('show');
          } else if(rechoice_model.get('id') == 'jilhvc5p2cy0atls8659a1fggjvvkmahwuspy2kr') {
            //-- The Answer we want
            $('.popover').popover('dispose');

            $relation_el = $('#selected-choice');
            $relation_el.popover({
              container: 'body',
              html: true,
              animation: false,
              content: function() {
                return "Correct, click 'submit' at the bottom of this page.";
              },
              placement: 'top'
            });
            $relation_el.popover('show');

            $('#submit_button').show();
            $next_el = $('#submit_button');
            $next_el.popover('hide');
            $next_el.popover({
              container: 'body',
              html: true,
              animation: false,
              content: function() {
                return 'Remember, you can always go back to your previous selection by clicking on the text of your current selection. <br />Click on \'submit\' to continue'
              },
              placement: 'top'
            });
            $next_el.popover('show');
          } else {
            this.better_answer();
          }
          this.event_counter++;
        },

        onDestroy: function() {
          $('.popover').popover('dispose');
        },

        onStartUpLogic: function() {
          console.log('onStartUpLogic: 1 0 1');

          $relates_el = $('ul a.list-group-item:nth(0)');
          $relates_el.popover({
            container: 'body',
            html: true,
            animation: false,
            content: function() {
              return 'Determine the relationship based ONLY on the text included. Click on the \'relates to\' box';
            },
            placement: 'top'
          });
          $relates_el.popover('show');
        }

      }
    },
  },
  3: {
    1: {
      incorrectClick: function() {
        this.better_answer();
      },
      choiceClick: function(rechoice_model) {
        if(rechoice_model.get('id') == '04110gzdcxz8niuv83ev08ut7lv0xep4iym5sxm5') {
          this.correct_answer();
        } else if(rechoice_model.get('id') == 'qq84lkjfh46gmx4a9n1jpwxwrmbajsy1qctb9u8j') {
          $('.popover').popover('dispose');
          this.hide_submit()
        } else {
          this.better_answer();
        }
      }
    },

    2: {
      incorrectClick: function() {
        this.better_answer();
      },
      choiceClick: function(rechoice_model) {
        if(rechoice_model.get('id') == 'u0q779rcrevnu6aki694dqka4fnfwvwqgpl06ybl') {
          this.correct_answer();
        } else if(rechoice_model.get('id') == 'qq84lkjfh46gmx4a9n1jpwxwrmbajsy1qctb9u8j') {
          $('.popover').popover('dispose');
          this.hide_submit()
        } else {
          this.better_answer();
        }
      }
    },

    3: {
      incorrectClick: function() {
        this.better_answer();
      },
      choiceClick: function(rechoice_model) {
        if(rechoice_model.get('id') == 'rhkmksv5jh0vn7p47uk3fwdior6mlgaubwh1l6ow') {
          this.correct_answer();
        } else if(rechoice_model.get('id') == 'qq84lkjfh46gmx4a9n1jpwxwrmbajsy1qctb9u8j') {
          $('.popover').popover('dispose');
          this.hide_submit()
        } else {
          this.better_answer();
        }
      }
    }
  },

  4: {
    1: {
      incorrectClick: function() {
        this.better_answer();
      },
      choiceClick: function(rechoice_model) {
        if(rechoice_model.get('id') == 'yrrb92b8vtjmcagjj4nx43sbj8wey2moqagk9ea5') {
          this.correct_answer();
        } else if(rechoice_model.get('id') == '8qota4u8hwtcyp65kz9zm0vjyuxwjt12sko084sn') {
          $('.popover').popover('dispose');
          this.hide_submit()
        } else {
          this.better_answer();
        }
      }
    },

    2: {
      incorrectClick: function() {
        this.better_answer();
      },
      choiceClick: function(rechoice_model) {
        if(rechoice_model.get('id') == 'lt18qfxd1ehj7ymxb29wrv6qa41mocwe6eor9dna') {
          this.correct_answer();
        } else if(rechoice_model.get('id') == '8qota4u8hwtcyp65kz9zm0vjyuxwjt12sko084sn') {
          $('.popover').popover('dispose');
          this.hide_submit()
        } else {
          this.better_answer();
        }
      }
    },

    3: {
      incorrectClick: function() {
        this.better_answer();
      },
      choiceClick: function(rechoice_model) {
        if(rechoice_model.get('id') == 'jyiczzhhupcp7cmebl422ax5dxe1jkwuq647oaw2') {
          this.correct_answer();
        } else if(rechoice_model.get('id') == '8qota4u8hwtcyp65kz9zm0vjyuxwjt12sko084sn') {
          $('.popover').popover('dispose');
          this.hide_submit()
        } else {
          this.better_answer();
        }
      }
    }
  },

  5: {
    1: {
      incorrectClick: function() {
        this.better_answer();
      },
      choiceClick: function(rechoice_model) {
        if(rechoice_model.get('id') == 'xdnolju6wacvakqnmz6237zwbh0ta3ftw8mdcp50') {
          this.correct_answer();
        } else if(rechoice_model.get('id') == 'txh8mu2mrnrffik893gr5h0ir7b1y7plgw94n4j7') {
          $('.popover').popover('dispose');
          this.hide_submit()
        } else {
          this.better_answer();
        }
      }
    },

    2: {
      incorrectClick: function() {
        this.better_answer();
      },
      choiceClick: function(rechoice_model) {
        if(rechoice_model.get('id') == 'am1wc2yvdcvwcb3yi298xesplbdktzku6wis49iw') {
          this.correct_answer();
        } else if(rechoice_model.get('id') == 'txh8mu2mrnrffik893gr5h0ir7b1y7plgw94n4j7') {
          $('.popover').popover('dispose');
          this.hide_submit()
        } else {
          this.better_answer();
        }
      }
    },

    3: {
      incorrectClick: function() {
        this.better_answer();
      },
      choiceClick: function(rechoice_model) {
        if(rechoice_model.get('id') == '5ex6vuro19zeneiwlc8yze6dsq1coxvlpojolwgy') {
          this.correct_answer();
        } else if(rechoice_model.get('id') == 'txh8mu2mrnrffik893gr5h0ir7b1y7plgw94n4j7') {
          $('.popover').popover('dispose');
          this.hide_submit()
        } else {
          this.better_answer();
        }
      }
    },

    4: {
      incorrectClick: function() {
        this.better_answer();
      },
      choiceClick: function(rechoice_model) {
        if(rechoice_model.get('id') == 'xdnolju6wacvakqnmz6237zwbh0ta3ftw8mdcp50') {
          this.correct_answer();
        } else {
          this.better_answer();
        }
      }
    },

  }

});


RETrainingAction = Backbone.Marionette.View.extend({
  /* Where all the custom UI testing is going down
   * Intentionally does not have a model or collection because option.training_data
   * and option.training_rule may come from an Step or an Instruction
   *
   * this.model = None
   * this.collection = None
   */
  template: '#training-action-template',

  initialize: function() {
    this.event_counter = 0;

    try {
      var pos = _.reject(this.getOption('position'), _.isNull);
      if(pos.length == 2) {
        var action_version = RETrainingActionLogic[pos[0]||0][pos[1]||0] || {};
      } else if (pos.length == 3) {
        var action_version = RETrainingActionLogic[pos[0]||0][pos[1]||0][pos[2]] || {};
      }
      _.extend(this, action_version);
    }
    catch(err) {};

    this.listenTo(channel, 'training:re:selectedconcept:clear', function() { this.clearClick(); });
    this.listenTo(channel, 'training:re:concept:incorrect:click', function(reconcept) { this.incorrectClick(reconcept); });
    this.listenTo(channel, 'training:re:choice:click', function(rechoice_model) { this.choiceClick(rechoice_model); });

    this.listenTo(channel, 'training:next:instruction', function() { this.onStartUpLogic(); })
    this.listenTo(channel, 'training:set:instruction', function() { this.onStartUpLogic(); })
  },

  onAttach: function() {
    var self = this;

    if(this.getOption('training_data')) {
      var RETreeApp = Backbone.Marionette.Application.extend({
        region: '#tree-action-area',

        onStart: function() {
          var main = this.getRegion();
          main.show( new Tree({
            'mode': 're',
            'training_data': self.getOption('training_data')
          }));
        }
      });
      var tree_app = new RETreeApp();
      tree_app.start();
      self.onStartUpLogic();
    }

  },

  onDestroy: function() {
    $('.popover').popover('dispose');
  },

  //----------------
  onStartUpLogic: function() {},
  clearClick: function() {
    $('.popover').popover('dispose');
    this.hide_submit();
  },
  incorrectClick: function(reconcept) {},
  choiceClick: function(rechoice_model) {},

  better_answer: function() {
    $('.popover').popover('dispose');
    this.hide_submit();
    $relation_el = $('#selected-choice');
    $relation_el.popover({
      container: 'body',
      html: true,
      animation: false,
      content: function() {
        return 'There is a better answer. Try again!';
      },
      placement: 'top'
    });
    $relation_el.popover('show');
  },

  show_submit: function() {
    $('#submit_button').removeClass('disabled');
    $('#submit_button').prop('disabled', false);
  },

  hide_submit: function() {
    $('#submit_button').addClass('disabled');
    $('#submit_button').prop('disabled', true);
  },

  correct_answer: function() {
    $('.popover').popover('dispose');
    this.show_submit();
    $relation_el = $('#selected-choice');
    $relation_el.popover({
      container: 'body',
      html: true,
      animation: false,
      content: function() {
        return "Correct, click 'submit' at the bottom of this page.";
      },
      placement: 'top'
    });
    $relation_el.popover('show');
  }

});


RETrainingStepView = TrainingStepView.extend({

  onRender: function() {
    this.showChildView('text', new TrainingStepTextView({'model': this.model}));
    this.showChildView('footer', new TrainingFooterView({'model': this.model}));
  }

});


RETrainingModuleView = TrainingModuleView.extend({

  onRender: function() {
    if(this.model.get('steps').length > 1) {
      this.showChildView('progress', new TrainingStepProgress({'collection': this.model.get('steps')}));
    }

    var step = this.model.get('steps').get_active();
    this.showChildView('step', new RETrainingStepView({'model': step}));
  }

});


RETrainingTaskView = TrainingTaskView.extend({

  initialize: function() {
    this.collection = new TrainingModuleCollection([], {'task_type': 're'});
    this.collection.fetch();
  },

  onRender: function() {
    if(this.collection.length) {
      this.showChildView('module-navigation', new TrainingModuleNavigation({'collection': this.collection}));

      var module = this.collection.get_active();
      this.showChildView('module', new RETrainingModuleView({'model': module}));

    }
  }

});

