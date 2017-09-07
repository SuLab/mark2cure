var channel = Backbone.Radio.channel('mark2cure');

//
// Views
//



RETrainingActionLogic = Object({

  //-- Modules
  1: {},
  2: {},
  3: {
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

        onStartUpLogic: function() {
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
    1: {}
  },
  4: {}
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

    var pos = this.getOption('position');

    var action_version = RETrainingActionLogic[pos[0]][pos[1]][pos[2]];
    _.extend(this, action_version);

    this.listenTo(channel, 'training:re:selectedconcept:clear', function() { this.clearClick(); });
    this.listenTo(channel, 'training:re:concept:incorrect:click', function(reconcept) { this.incorrectClick(reconcept); });
    this.listenTo(channel, 'training:re:choice:click', function(rechoice_model) { this.choiceClick(rechoice_model); });

  },


  //----------------
  onStartUpLogic: function() {},
  clearClick: function() {},
  incorrectClick: function(reconcept) {},
  choiceClick: function(rechoice_model) {},

  better_answer: function() {
    $('.popover').popover('dispose');
    $('#submit_button').hide();
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

});



RETrainingStepView = TrainingStepView.extend({

  onRender: function() {
    console.log('RETrainingStepView', this);
    this.showChildView('text', new TrainingStepTextView({'model': this.model}));
    this.showChildView('footer', new TrainingFooterView({'model': this.model}));
  }

});


RETrainingModuleView = TrainingModuleView.extend({

  onRender: function() {
    this.showChildView('progress', new TrainingStepProgress({'collection': this.model.get('steps')}));
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

