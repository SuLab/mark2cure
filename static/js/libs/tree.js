/*
 *  Models & Collections
 */
Relation = Backbone.RelationalModel.extend({
  defaults: {
    id: '',
    text: '',
    selected: false
  },

  relations: [{
    type: 'HasMany',
    key: 'children',

    relatedModel: 'Relation',
    collectionType: 'RelationList',

    reverseRelation : {
      key : 'parentRelation',
      includeInJSON: false,
    }
  }],

  get_selected: function() {
    /* Try to find an instance of 'this' model type in the store */
    var model = Backbone.Relational.store.find( this, {"selected": true});
    // console.log('model:', model);

    if ( !model && _.isObject( attributes ) ) {
      var coll = Backbone.Relational.store.getCollection( this );

      model = coll.find( function( m ) {
        return m.selected === true;
      });
    }
    return model;
  }

});

RelationList = Backbone.Collection.extend({
  model: Relation,
  url: '/api/v1/words',
});

function get_stype(the_current_relation) {
  return $(this).hasClass('.gene');
};
/*
 * Views
 */
RelationView = Backbone.Marionette.ItemView.extend({
  template: _.template('<%= text %>'),
  tagName: 'a',
  className: 'list-group-item',

  events : {
    'mousedown' : 'mousedown',
  },

  mousedown : function(evt) {
    /*
     * 1. Set the current choice to the ID
     * 2. Set the backbutton reference if available
     * 3. Set the list updated to any children
     *
     */
    var children = this.model.get('children');
    Tree['convoChannel'].trigger('click', {'collection': children, 'choice': this.model});
  },
});


RelationCompositeView = Backbone.Marionette.CompositeView.extend({
  template: '#tree-template',
  templateHelpers: function() {
    return this.options.concepts
  },

  childView  : RelationView,
  childViewContainer: "ul",
  /*
  tagName   : 'div',
  className : 'paragraph',
  */

  ui: {
    'c1': '#c1',
    'c2': '#c2',
    'relation': '#relation',
    'c1_not_correct': '#c1_not_correct',
    'c2_not_correct': '#c2_not_correct',
    'incorrect_buttons': '.fa.fa-times-circle',
    'submit_button': '#submit_button'
  },

  events : {
    'mouseover @ui.incorrect_buttons': 'hoverIncorrect',
    'mouseout @ui.incorrect_buttons': 'hoverOutIncorrect',
    'mousedown @ui.relation': 'resetRelationship',
    'mousedown @ui.c1_not_correct': 'c1Error',
    'mousedown @ui.c2_not_correct': 'c2Error',
    'mousedown @ui.submit_button': 'submit_fade_in',
  },

  hoverIncorrect: function(evt) {
    $(evt.target).parent().addClass('incorrect');
  },
  hoverOutIncorrect: function(evt) {
    $(evt.target).parent().removeClass('incorrect');
  },

  resetRelationship: function(evt) {
    Tree['convoChannel'].trigger('back', this.options);
  },

  c1Error: function(evt) {
    Tree['convoChannel'].trigger('error', 'c_1');
  },
  c2Error: function(evt) {
    Tree['convoChannel'].trigger('error', 'c_2');
  },

  onRender : function() {
    var self = this;
    var choice = this.options['choice'];
    var concepts = this.options['concepts'];
    /*
     * console.log('[RelationCompositeView onRender] Choice:', choice);
     */

    if(this.options.first_draw) {
      if( _.has(concepts['c1'], 'fadeIn') ) {
        this.ui.c1.addClass('fade-in one');
      };
      if( _.has(concepts['c2'], 'fadeIn') ) {
        this.ui.c2.addClass('fade-in one');
      };
    }

    if(choice) {
      function get_stype_word(stype){
        if (stype == 'g'){
          stype_word = 'gene';
        } else if ( stype == 'c'){
          stype_word = 'drug';
        } else {
          stype_word = 'disease';
        }
        return stype_word;
      };

      if (choice.id == "zl4RlTGwZM9Ud3CCXpU2VZa7eQVnJj0MdbsRBMGy") {
        this.ui.relation.removeClass('disabled').text( concepts['c1'].text + choice.get('text') + get_stype_word(concepts['c1'].type) + ' concept' );

      } else if (choice.id == "RdKIrcaEOnM4DRk25g5jAfeNC6HSpsFZaiIPqZer") {
        this.ui.relation.removeClass('disabled').text( concepts['c2'].text + choice.get('text') + get_stype_word(concepts['c2'].type) + ' concept' );

        /* Training specific
         * (TODO) get out of here
         */
        if(concepts['c2']['text'] == 'Astrology') {
          this.ui.relation.removeClass('disabled').text( concepts['c2'].text + choice.get('text') + 'field of science');
        };

      } else {
        this.ui.relation.removeClass('disabled').text( choice.get('text') );
      };
    };

    if(choice || this.collection.parentRelation) {
      this.ui.relation.removeClass('disabled');
      self.ui.relation.addClass('relation-go-back');
    }

  }
});

Tree = new Backbone.Marionette.Application();
